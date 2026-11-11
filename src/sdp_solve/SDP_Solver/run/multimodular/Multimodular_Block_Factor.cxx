#include "Multimodular_Block_Factor.hxx"

#include "sdpb_util/Shared_Window_Array.hxx"
#include "sdpb_util/assert.hxx"

#include <algorithm>
#include <cstdint>
#include <vector>

namespace
{
  // count BigFloat values in an MPI-3 shared-memory window
  // (Shared_Window_Array), one record of 2 + limbs mp_limb_t each:
  // _mp_size, _mp_exp, then the limbs, as in an mpf_t of precision
  // El::gmp::Precision(). A rank writes a value through an mpf struct whose
  // limb pointer points into the window (view()), then stores its size and
  // exponent with publish(). After Fence(), other ranks read it with view().
  // The window may map at different addresses in different processes, so
  // each rank makes its own structs.
  class Shared_BigFloat_Array
  {
    Shared_Window_Array<mp_limb_t> window;
    int limbs; // mpf writes at most _mp_prec + 1 limbs
    size_t record;

    static int value_limbs()
    {
      mpf_t probe;
      mpf_init2(probe, El::gmp::Precision());
      const int result = probe->_mp_prec + 1;
      mpf_clear(probe);
      return result;
    }

  public:
    Shared_BigFloat_Array(const El::mpi::Comm &comm, const size_t count)
        : window(comm, count * (2 + value_limbs())),
          limbs(value_limbs()),
          record(2 + limbs)
    {
      static_assert(sizeof(mp_limb_t) == sizeof(long)
                      && sizeof(mp_limb_t) >= sizeof(int),
                    "size and exponent are stored as limbs");
      if(El::mpi::Rank(comm) == 0)
        for(size_t index = 0; index < count; ++index)
          {
            window[index * record] = 0;
            window[index * record + 1] = 0;
          }
      Fence();
    }

    // An mpf struct for the value at index: limbs in the window, size and
    // exponent as last published
    __mpf_struct view(const size_t index)
    {
      __mpf_struct v;
      v._mp_prec = limbs - 1;
      v._mp_size = static_cast<int>(static_cast<long>(window[index * record]));
      v._mp_exp = static_cast<long>(window[index * record + 1]);
      v._mp_d = &window[index * record + 2];
      return v;
    }
    // Store size and exponent of v, whose limbs are those of view(index)
    void publish(const size_t index, const __mpf_struct &v)
    {
      ASSERT(v._mp_d == &window[index * record + 2]);
      window[index * record] = static_cast<mp_limb_t>(long(v._mp_size));
      window[index * record + 1] = static_cast<mp_limb_t>(long(v._mp_exp));
    }
    void Fence() const { window.Fence(); }
  };

  // Lower triangle of an n x n BigFloat matrix in a Shared_BigFloat_Array
  // on the block's grid (all its ranks on one node): one copy per node
  // instead of one per rank. Each rank has its own El::Matrix of mpf structs
  // whose limb pointers point into the window. A rank writes elements
  // through its structs, then publish() copies their size and exponent to
  // the window; after sync(), pull() refreshes the other ranks' structs.
  // Strictly upper elements are never used: they point to a dummy zero.
  class Shared_Lower_Matrix
  {
    const int n;
    const int rank;
    Shared_BigFloat_Array array;
    El::BigFloat *structs = nullptr;
    mp_limb_t dummy[2] = {0, 0};

    size_t index(const int i, const int j) const
    { return size_t(j) * n - size_t(j) * (j - 1) / 2 + (i - j); }
    __mpf_struct *element(const int i, const int j)
    { return structs[i + size_t(j) * n].gmp_float.get_mpf_t(); }

  public:
    El::Matrix<El::BigFloat> matrix;

    Shared_Lower_Matrix(const int n, const El::mpi::Comm &comm)
        : n(n), rank(El::mpi::Rank(comm)), array(comm, size_t(n) * (n + 1) / 2)
    {
      static_assert(sizeof(El::BigFloat) == sizeof(__mpf_struct),
                    "El::BigFloat must be a bare mpf");
      // mpf structs filled in place (no constructors: limbs live in the
      // window, so ~BigFloat must never run on them)
      structs = static_cast<El::BigFloat *>(
        ::operator new(sizeof(El::BigFloat) * size_t(n) * n));
      for(int j = 0; j < n; ++j)
        for(int i = 0; i < n; ++i)
          {
            __mpf_struct *s = element(i, j);
            if(i >= j)
              *s = array.view(index(i, j));
            else
              {
                s->_mp_prec = 1;
                s->_mp_size = 0;
                s->_mp_exp = 0;
                s->_mp_d = dummy;
              }
          }
      matrix.Attach(n, n, structs, n);
    }
    ~Shared_Lower_Matrix()
    {
      matrix.Empty();
      ::operator delete(structs);
    }
    Shared_Lower_Matrix(const Shared_Lower_Matrix &) = delete;
    Shared_Lower_Matrix &operator=(const Shared_Lower_Matrix &) = delete;

    int Rank() const { return rank; }

    void publish(const int i, const int j)
    { array.publish(index(i, j), *element(i, j)); }
    void pull(const int i, const int j)
    {
      const __mpf_struct v = array.view(index(i, j));
      __mpf_struct *s = element(i, j);
      s->_mp_size = v._mp_size;
      s->_mp_exp = v._mp_exp;
    }
    // publish/pull the lower part of rows [i0, i1) x columns [j0, j1)
    void publish(const int i0, const int i1, const int j0, const int j1)
    {
      for(int j = j0; j < j1; ++j)
        for(int i = std::max(i0, j); i < i1; ++i)
          publish(i, j);
    }
    void pull(const int i0, const int i1, const int j0, const int j1)
    {
      for(int j = j0; j < j1; ++j)
        for(int i = std::max(i0, j); i < i1; ++i)
          pull(i, j);
    }
    void sync() { array.Fence(); }
  };

  // n x n El::Matrix<BigFloat> whose strictly upper elements have minimal
  // precision: they are never used, and this halves the memory of the
  // replicated factor (El::Matrix::Resize would allocate them in full).
  class Lower_BigFloat_Matrix
  {
    El::BigFloat *data = nullptr;
    size_t size = 0;

  public:
    El::Matrix<El::BigFloat> matrix;
    explicit Lower_BigFloat_Matrix(const int n) : size(size_t(n) * n)
    {
      data = static_cast<El::BigFloat *>(
        ::operator new(sizeof(El::BigFloat) * size));
      for(int j = 0; j < n; ++j)
        for(int i = 0; i < n; ++i)
          new(data + i + size_t(j) * n)
            El::BigFloat(0, i >= j ? El::gmp::Precision() : 1);
      matrix.Attach(n, n, data, n);
    }
    ~Lower_BigFloat_Matrix()
    {
      matrix.Empty();
      for(size_t k = 0; k < size; ++k)
        data[k].~BigFloat();
      ::operator delete(data);
    }
    Lower_BigFloat_Matrix(const Lower_BigFloat_Matrix &) = delete;
    Lower_BigFloat_Matrix &operator=(const Lower_BigFloat_Matrix &) = delete;
  };

  // True if all ranks of comm share a node (MPI shared memory works)
  bool on_one_node(const El::mpi::Comm &comm)
  {
    MPI_Comm node;
    MPI_Comm_split_type(comm.comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL,
                        &node);
    int node_size, size;
    MPI_Comm_size(node, &node_size);
    MPI_Comm_size(comm.comm, &size);
    MPI_Comm_free(&node);
    return node_size == size;
  }

  // Allgather a [STAR,VR] matrix M in 64-column panels (small transient
  // buffers) and call f(i, j, M(i,j)) on every rank; if lower, only i >= j.
  template <class F>
  void
  allgather_columns(const El::DistMatrix<El::BigFloat, El::STAR, El::VR> &M,
                    const bool lower, const F &f)
  {
    const int height = M.Height(), width = M.Width();
    const int panel = 64;
    for(int j0 = 0; j0 < width; j0 += panel)
      {
        const int j1 = std::min(width, j0 + panel);
        const int r0 = lower ? std::min(j0, height) : 0;
        El::DistMatrix<El::BigFloat, El::STAR, El::STAR> part(
          M(El::IR(r0, height), El::IR(j0, j1)));
        for(int j = j0; j < j1; ++j)
          for(int i = lower ? j : 0; i < height; ++i)
            f(i, j, part.GetLocalCRef(i - r0, j - j0));
      }
  }

  // Cholesky (lower, in place) of the n x n matrix A, a view of a matrix
  // that every rank of the grid holds in full, with the same values.
  // In each recursion level, A = [A11 *; A21 A22]:
  //   A11 := chol(A11)                       (recursive)
  //   A21 := L21, L21^T = L11^{-1} A21^T     (columns of A21^T split)
  //   A22 -= L21 L21^T, lower triangle       (columns of A22 split)
  //   A22 := chol(A22)                       (recursive)
  // The split parts are computed by each rank for its columns
  // ([STAR,VR]: column j on rank j mod g) and then allgathered.
  // Below config.distributed_min_size (or on one rank) every rank computes
  // the whole factorization.
  void
  distributed_multimodular_cholesky(El::Matrix<El::BigFloat> &A,
                                    const El::Grid &grid,
                                    const Multimodular_Blas_Config &config)
  {
    const int n = A.Height();
    Multimodular_Blas_Timings timings;
    if(grid.Size() == 1 || n < config.distributed_min_size)
      {
        multimodular_cholesky_lower(A, config, timings);
        return;
      }
    const int n1 = n / 2, n2 = n - n1;
    const auto top = El::IR(0, n1), bottom = El::IR(n1, n);
    auto A11 = A(top, top);
    auto A21 = A(bottom, top);
    auto A22 = A(bottom, bottom);

    distributed_multimodular_cholesky(A11, grid, config);

    // L21^T = L11^{-1} A21^T
    {
      El::DistMatrix<El::BigFloat, El::STAR, El::VR> L21_T(n1, n2, grid);
      auto &local = L21_T.Matrix();
      for(int jLoc = 0; jLoc < L21_T.LocalWidth(); ++jLoc)
        {
          const int j = L21_T.GlobalCol(jLoc);
          for(int i = 0; i < n1; ++i)
            local.Set(i, jLoc, A21.CRef(j, i));
        }
      multimodular_trsm_left_lower(A11, local, config, timings);
      allgather_columns(L21_T, false,
                        [&](int i, int j, const El::BigFloat &value) {
                          A21.Set(j, i, value);
                        });
    }

    // A22 -= L21 L21^T, lower triangle
    {
      El::DistMatrix<El::BigFloat, El::STAR, El::VR> A22_columns(n2, n2, grid);
      auto &local = A22_columns.Matrix();
      const int width = A22_columns.LocalWidth();
      El::Matrix<El::BigFloat> L21_rows(width, n1);
      std::vector<int> row_start(width);
      for(int jLoc = 0; jLoc < width; ++jLoc)
        {
          const int j = A22_columns.GlobalCol(jLoc);
          row_start[jLoc] = j;
          for(int i = j; i < n2; ++i)
            local.Set(i, jLoc, A22.CRef(i, j));
          for(int l = 0; l < n1; ++l)
            L21_rows.Set(jLoc, l, A21.CRef(j, l));
        }
      multimodular_gemm_sub(local, A21, El::TRANSPOSE, L21_rows, row_start,
                            config, timings);
      allgather_columns(A22_columns, true,
                        [&](int i, int j, const El::BigFloat &value) {
                          A22.Set(i, j, value);
                        });
    }

    distributed_multimodular_cholesky(A22, grid, config);
  }

  // S := lower triangle of the distributed `schur`: each rank writes the
  // entries it owns.
  void fill_shared(const El::DistMatrix<El::BigFloat> &schur,
                   Shared_Lower_Matrix &S)
  {
    const int n = schur.Height();
    for(int jLoc = 0; jLoc < schur.LocalWidth(); ++jLoc)
      for(int iLoc = 0; iLoc < schur.LocalHeight(); ++iLoc)
        {
          const int i = schur.GlobalRow(iLoc), j = schur.GlobalCol(jLoc);
          if(i < j)
            continue;
          mpf_set(S.matrix.Ref(i, j).gmp_float.get_mpf_t(),
                  schur.GetLocalCRef(iLoc, jLoc).gmp_float.get_mpf_t());
          S.publish(i, j);
        }
    S.sync();
    S.pull(0, n, 0, n);
  }

  // Cholesky of rows/columns [r0, r0+n) of the shared matrix, as
  // distributed_multimodular_cholesky but writing in place instead of
  // allgathering. The parts below config.distributed_min_size are computed by
  // rank 0 of the grid only.
  void shared_multimodular_cholesky(Shared_Lower_Matrix &S, const int r0,
                                    const int n, const int grid_size,
                                    const Multimodular_Blas_Config &config)
  {
    Multimodular_Blas_Timings timings;
    const int rank = S.Rank();
    if(grid_size == 1 || n < config.distributed_min_size)
      {
        if(rank == 0)
          {
            auto A = S.matrix(El::IR(r0, r0 + n), El::IR(r0, r0 + n));
            multimodular_cholesky_lower(A, config, timings);
            S.publish(r0, r0 + n, r0, r0 + n);
          }
        S.sync();
        if(rank != 0)
          S.pull(r0, r0 + n, r0, r0 + n);
        return;
      }
    const int n1 = n / 2, n2 = n - n1;
    shared_multimodular_cholesky(S, r0, n1, grid_size, config);
    const auto rows1 = El::IR(r0, r0 + n1), rows2 = El::IR(r0 + n1, r0 + n);
    const auto A11 = S.matrix(rows1, rows1);
    auto A21 = S.matrix(rows2, rows1);
    auto A22 = S.matrix(rows2, rows2);

    // This rank's rows j of A21 (j = rank mod grid_size)
    std::vector<int> mine;
    for(int j = rank; j < n2; j += grid_size)
      mine.push_back(j);
    const int w = mine.size();

    // L21^T = L11^{-1} A21^T for this rank's rows, written in place
    {
      El::Matrix<El::BigFloat> T(n1, w);
      for(int jj = 0; jj < w; ++jj)
        for(int i = 0; i < n1; ++i)
          T.Set(i, jj, A21.CRef(mine[jj], i));
      multimodular_trsm_left_lower(A11, T, config, timings);
      for(int jj = 0; jj < w; ++jj)
        for(int i = 0; i < n1; ++i)
          {
            A21.Set(mine[jj], i, T.CRef(i, jj));
            S.publish(r0 + n1 + mine[jj], r0 + i);
          }
      S.sync();
      S.pull(r0 + n1, r0 + n, r0, r0 + n1);
    }

    // A22 -= L21 L21^T, this rank's columns, lower part, written in place
    {
      El::Matrix<El::BigFloat> C(n2, w), L21_rows(w, n1);
      std::vector<int> row_start(w);
      for(int jj = 0; jj < w; ++jj)
        {
          const int j = mine[jj];
          row_start[jj] = j;
          for(int i = j; i < n2; ++i)
            C.Set(i, jj, A22.CRef(i, j));
          for(int l = 0; l < n1; ++l)
            L21_rows.Set(jj, l, A21.CRef(j, l));
        }
      multimodular_gemm_sub(C, A21, El::TRANSPOSE, L21_rows, row_start, config,
                            timings);
      for(int jj = 0; jj < w; ++jj)
        {
          const int j = mine[jj];
          for(int i = j; i < n2; ++i)
            {
              A22.Set(i, j, C.CRef(i, jj));
              S.publish(r0 + n1 + i, r0 + n1 + j);
            }
        }
      S.sync();
      S.pull(r0 + n1, r0 + n, r0 + n1, r0 + n);
    }

    shared_multimodular_cholesky(S, r0 + n1, n2, grid_size, config);
  }

  // L := the lower triangle of the distributed matrix `schur`, on every rank
  // of its grid. Gathered in column panels, rows j0..n of columns j0..j1 at
  // a time, so the transient buffers stay small and only the lower triangle
  // moves.
  void gather_lower(const El::DistMatrix<El::BigFloat> &schur,
                    El::Matrix<El::BigFloat> &L)
  {
    const int n = schur.Height();
    ASSERT_EQUAL(L.Height(), n);
    const int panel = 64;
    for(int j0 = 0; j0 < n; j0 += panel)
      {
        const int j1 = std::min(n, j0 + panel);
        El::DistMatrix<El::BigFloat, El::STAR, El::STAR> part(
          schur(El::IR(j0, n), El::IR(j0, j1)));
        for(int j = j0; j < j1; ++j)
          for(int i = j; i < n; ++i)
            mpf_set(L.Ref(i, j).gmp_float.get_mpf_t(),
                    part.GetLocalCRef(i - j0, j - j0).gmp_float.get_mpf_t());
      }
  }

  // cholesky := schur, with its lower triangle replaced by L's
  void set_local(const El::DistMatrix<El::BigFloat> &schur,
                 const El::Matrix<El::BigFloat> &L,
                 El::DistMatrix<El::BigFloat> &cholesky)
  {
    cholesky = schur;
    for(int jLoc = 0; jLoc < cholesky.LocalWidth(); ++jLoc)
      for(int iLoc = 0; iLoc < cholesky.LocalHeight(); ++iLoc)
        {
          const int i = cholesky.GlobalRow(iLoc), j = cholesky.GlobalCol(jLoc);
          if(i >= j)
            cholesky.SetLocal(iLoc, jLoc, L.CRef(i, j));
        }
  }

  // B := L^{-1} B for a grid on one node, through a shared-memory window
  // holding B (one copy per block): each rank writes the entries of B it
  // owns ([MC,MR]), solves its columns (j = rank mod g) in place in the
  // window, and reads back the entries it owns. No [STAR,VR] copy and no El
  // redistribution buffers.
  void solve_shared(const El::Matrix<El::BigFloat> &L,
                    El::DistMatrix<El::BigFloat> &B,
                    const Multimodular_Blas_Config &config, Timers &timers)
  {
    Multimodular_Blas_Timings timings;
    const int P = B.Height(), N = B.Width();
    const El::mpi::Comm &comm = B.Grid().Comm();
    const int rank = El::mpi::Rank(comm), grid_size = El::mpi::Size(comm);

    Scoped_Timer to_shared_timer(timers, "to_shared");
    Shared_BigFloat_Array window(comm, size_t(P) * N);
    for(int jLoc = 0; jLoc < B.LocalWidth(); ++jLoc)
      for(int iLoc = 0; iLoc < B.LocalHeight(); ++iLoc)
        {
          const size_t index
            = B.GlobalRow(iLoc) + size_t(B.GlobalCol(jLoc)) * P;
          __mpf_struct v = window.view(index);
          mpf_set(&v, B.GetLocalCRef(iLoc, jLoc).gmp_float.get_mpf_t());
          window.publish(index, v);
        }
    window.Fence();
    to_shared_timer.stop();

    {
      Scoped_Timer timer(timers, "multimodular");
      std::vector<int> columns;
      for(int j = rank; j < N; j += grid_size)
        columns.push_back(j);
      const int w = columns.size();
      // mpf structs of this rank's columns, pointing into the window
      El::BigFloat *structs = static_cast<El::BigFloat *>(::operator new(
        sizeof(El::BigFloat) * std::max<size_t>(1, size_t(P) * w)));
      for(int jj = 0; jj < w; ++jj)
        for(int i = 0; i < P; ++i)
          *structs[i + size_t(jj) * P].gmp_float.get_mpf_t()
            = window.view(i + size_t(columns[jj]) * P);
      {
        El::Matrix<El::BigFloat> X;
        X.Attach(P, w, structs, P);
        multimodular_trsm_left_lower(L, X, config, timings);
      }
      for(int jj = 0; jj < w; ++jj)
        for(int i = 0; i < P; ++i)
          window.publish(i + size_t(columns[jj]) * P,
                         *structs[i + size_t(jj) * P].gmp_float.get_mpf_t());
      ::operator delete(structs);
      window.Fence();
    }

    Scoped_Timer from_shared_timer(timers, "from_shared");
    for(int jLoc = 0; jLoc < B.LocalWidth(); ++jLoc)
      for(int iLoc = 0; iLoc < B.LocalHeight(); ++iLoc)
        {
          const __mpf_struct v
            = window.view(B.GlobalRow(iLoc) + size_t(B.GlobalCol(jLoc)) * P);
          mpf_set(B.Matrix().Ref(iLoc, jLoc).gmp_float.get_mpf_t(), &v);
        }
  }

  // B := L^{-1} B. The columns of B are split among the ranks of the grid
  // ([STAR,VR] distribution); each rank solves for its own columns with
  // its full copy of L.
  void solve_columns(const El::Matrix<El::BigFloat> &L,
                     El::DistMatrix<El::BigFloat> &B,
                     const Multimodular_Blas_Config &config, Timers &timers)
  {
    Scoped_Timer to_columns_timer(timers, "to_columns");
    El::DistMatrix<El::BigFloat, El::STAR, El::VR> B_columns(B);
    to_columns_timer.stop();
    {
      Scoped_Timer timer(timers, "multimodular");
      Multimodular_Blas_Timings timings;
      multimodular_trsm_left_lower(L, B_columns.Matrix(), config, timings);
    }
    Scoped_Timer timer(timers, "from_columns");
    B = B_columns;
  }
}

struct Multimodular_Block_Factor::Storage
{
  std::unique_ptr<Shared_Lower_Matrix> shared;
  std::unique_ptr<Lower_BigFloat_Matrix> local;
  El::Matrix<El::BigFloat> &L()
  { return shared ? shared->matrix : local->matrix; }
};

Multimodular_Block_Factor::Multimodular_Block_Factor(
  const El::DistMatrix<El::BigFloat> &schur,
  const Multimodular_Blas_Config &config, const bool allow_shared)
    : storage(std::make_unique<Storage>()), config(config)
{
  ASSERT_EQUAL(schur.Height(), schur.Width());
  const El::mpi::Comm &comm = schur.Grid().Comm();
  if(allow_shared && on_one_node(comm))
    storage->shared
      = std::make_unique<Shared_Lower_Matrix>(schur.Height(), comm);
  else
    storage->local = std::make_unique<Lower_BigFloat_Matrix>(schur.Height());
}

Multimodular_Block_Factor::~Multimodular_Block_Factor() = default;

bool Multimodular_Block_Factor::shared() const
{ return storage->shared != nullptr; }

void Multimodular_Block_Factor::cholesky(
  const El::DistMatrix<El::BigFloat> &schur,
  El::DistMatrix<El::BigFloat> &cholesky, Timers &timers)
{
  if(shared())
    {
      Shared_Lower_Matrix &S = *storage->shared;
      {
        Scoped_Timer timer(timers, "fill_shared");
        fill_shared(schur, S);
      }
      {
        Scoped_Timer timer(timers, "multimodular");
        shared_multimodular_cholesky(S, 0, schur.Height(), schur.Grid().Size(),
                                     config);
      }
    }
  else
    {
      El::Matrix<El::BigFloat> &L = storage->L();
      {
        Scoped_Timer timer(timers, "allgather");
        gather_lower(schur, L);
      }
      Scoped_Timer timer(timers, "multimodular");
      distributed_multimodular_cholesky(L, schur.Grid(), config);
    }
  Scoped_Timer timer(timers, "set_local");
  set_local(schur, storage->L(), cholesky);
}

void Multimodular_Block_Factor::solve(El::DistMatrix<El::BigFloat> &B,
                                      Timers &timers)
{
  const El::Matrix<El::BigFloat> &L = storage->L();
  ASSERT_EQUAL(B.Height(), L.Height());
  if(B.Grid().Size() == 1)
    {
      // In place
      Scoped_Timer timer(timers, "multimodular");
      Multimodular_Blas_Timings timings;
      multimodular_trsm_left_lower(L, B.Matrix(), config, timings);
    }
  else if(shared())
    solve_shared(L, B, config, timers);
  else
    solve_columns(L, B, config, timers);
}

size_t multimodular_block_factor_bytes(const int n, const int width,
                                       const int grid_size,
                                       const size_t bigfloat_bytes,
                                       const Multimodular_Blas_Config &config)
{
  if(n <= 0)
    return 0;
  // Primes for products of (bits + 32)-bit integers: at least 20 bits each
  // below the cap of Fmpz_Comb's primes (see blas_residues.hxx)
  const size_t bits = El::gmp::Precision() + config.guard_bits;
  const size_t num_primes = (2 * bits + config.growth_bits + 64) / 20 + 1;
  const size_t nn = size_t(n) * n, g = std::max(grid_size, 1);
  const size_t L_window = bigfloat_bytes * nn / 2;
  // Per rank: the structs (n x n), A21's residues (n/2 x n/2, every rank),
  // its C and L21 copies (n^2 / 2g), the residues of op(B) and C
  // (k + m = n rows, column_chunk columns)
  const size_t cholesky_rank
    = sizeof(El::BigFloat) * nn + sizeof(int32_t) * num_primes * nn / 4
      + bigfloat_bytes * nn / (2 * g)
      + sizeof(int32_t) * num_primes * n * config.column_chunk;
  // Shared copy of B, and per rank the X residues (at most
  // left_looking_bytes / 2) and the panel residues
  const size_t B_window = bigfloat_bytes * size_t(n) * width;
  const size_t solve_rank
    = std::min(sizeof(int32_t) * num_primes * n * width,
               config.left_looking_bytes / 2)
      + 2 * sizeof(int32_t) * num_primes * config.base_size * n;
  return L_window
         + std::max(g * cholesky_rank, B_window + g * solve_rank);
}
