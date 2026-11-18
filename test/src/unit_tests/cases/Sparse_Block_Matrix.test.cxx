#include <catch2/catch_amalgamated.hpp>
#include "sdp_solve/Block_Matrix/Sparse_Block_Matrix.hxx"
#include "unit_tests/util/util.hxx"

// Sparse_Block_Matrix keeps a DistMatrix's local nonzeros and gives it back
// with to_dense(). Matrices with zeros inside columns, whole zero columns and
// rows, on grids of 1 to all ranks, with several alignments: the round trip
// must give every element back exactly.
TEST_CASE("Sparse_Block_Matrix")
{
  const int world_size = El::mpi::Size(El::mpi::COMM_WORLD);
  const int world_rank = El::mpi::Rank(El::mpi::COMM_WORLD);
  const int height = GENERATE(1, 7, 30);
  const int width = GENERATE(1, 13);
  for(int grid_size = 1; grid_size <= world_size; ++grid_size)
    DYNAMIC_SECTION("height=" << height << " width=" << width
                              << " grid_size=" << grid_size)
    {
      // The same matrix on every rank: about 40% zeros, placed at random,
      // plus one zero column and one zero row
      El::InitializeRandom(true);
      El::Matrix<El::BigFloat> global
        = Test_Util::random_matrix(height, width);
      El::mpi::Broadcast(global.Buffer(), global.MemorySize(), 0,
                         El::mpi::COMM_WORLD);
      for(int j = 0; j < width; ++j)
        for(int i = 0; i < height; ++i)
          if((i * 7 + j * 3 + i * j) % 5 < 2 || j == width / 2
             || i == height / 3)
            global.Set(i, j, El::BigFloat(0));

      El::mpi::Comm comm;
      El::mpi::Split(El::mpi::COMM_WORLD, world_rank < grid_size ? 0 : 1,
                     world_rank, comm);
      if(world_rank < grid_size)
        {
          const El::Grid grid(comm);
          for(int align = 0; align < std::min(grid_size, 2); ++align)
            {
              CAPTURE(align);
              El::DistMatrix<El::BigFloat> dense(grid);
              dense.AlignCols(align % grid.Height(), false);
              dense.AlignRows(align % grid.Width(), false);
              dense.Resize(height, width);
              for(int jLoc = 0; jLoc < dense.LocalWidth(); ++jLoc)
                for(int iLoc = 0; iLoc < dense.LocalHeight(); ++iLoc)
                  dense.SetLocal(iLoc, jLoc,
                                 global.Get(dense.GlobalRow(iLoc),
                                            dense.GlobalCol(jLoc)));

              const Sparse_Block_Matrix sparse(dense);
              size_t nonzeros = 0;
              for(int jLoc = 0; jLoc < dense.LocalWidth(); ++jLoc)
                for(int iLoc = 0; iLoc < dense.LocalHeight(); ++iLoc)
                  nonzeros += dense.GetLocal(iLoc, jLoc) != El::BigFloat(0);
              REQUIRE(sparse.num_nonzeros_local() == nonzeros);

              El::DistMatrix<El::BigFloat> result(grid);
              sparse.to_dense(result);
              REQUIRE(result.Height() == height);
              REQUIRE(result.Width() == width);
              REQUIRE(result.ColAlign() == dense.ColAlign());
              REQUIRE(result.RowAlign() == dense.RowAlign());
              for(int jLoc = 0; jLoc < result.LocalWidth(); ++jLoc)
                for(int iLoc = 0; iLoc < result.LocalHeight(); ++iLoc)
                  {
                    CAPTURE(iLoc);
                    CAPTURE(jLoc);
                    REQUIRE(result.GetLocal(iLoc, jLoc)
                            == dense.GetLocal(iLoc, jLoc));
                  }
            }
        }
      El::mpi::Free(comm);
    }
}
