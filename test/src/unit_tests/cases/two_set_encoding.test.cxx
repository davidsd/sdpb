// pmp2sdp --encodingPoints (two-set encoding): the encoding-point rules and
// the second step (approximation points -> encoding points).

#include "pmp/Polynomial_Vector_Matrix.hxx"
#include "pmp2sdp/Dual_Constraint_Group.hxx"
#include "sdpb_util/Damped_Rational.hxx"
#include "sdpb_util/boost_serialization.hxx"

#include <catch2/catch_amalgamated.hpp>

#include "unit_tests/util/util.hxx"

#include <boost/math/constants/constants.hpp>

std::vector<Boost_Float>
sample_points(const size_t &num_points, const Damped_Rational &prefactor);

std::vector<Boost_Float>
sample_scalings(const std::vector<Boost_Float> &points,
                const Damped_Rational &damped_rational);

std::array<Polynomial_Vector, 2>
bilinear_basis(const std::vector<El::BigFloat> &sample_points,
               const std::vector<El::BigFloat> &sample_scalings);

std::array<El::Matrix<El::BigFloat>, 2>
sample_bilinear_basis(const std::array<Polynomial_Vector, 2> &bilinear_basis,
                      const std::vector<El::BigFloat> &sample_points,
                      const std::vector<El::BigFloat> &sample_scalings);

std::array<El::Matrix<El::BigFloat>, 2> sample_bilinear_basis_by_recurrence(
  const std::array<size_t, 2> &basis_sizes,
  const std::vector<El::BigFloat> &sample_points,
  const std::vector<El::BigFloat> &sample_scalings);

namespace
{
  Boost_Float rel_err(const Boost_Float &a, const Boost_Float &b)
  {
    return Boost_Float(abs(a - b)) / std::max(Boost_Float(abs(a)), Boost_Float(abs(b)));
  }

  // A deterministic "random" polynomial of the given degree
  Polynomial test_polynomial(const size_t degree, const int seed)
  {
    Polynomial p;
    p.coefficients.clear();
    for(size_t i = 0; i <= degree; ++i)
      {
        const int v = (37 * static_cast<int>(i) + 11 * seed + 5) % 23 - 11;
        p.coefficients.emplace_back(El::BigFloat(v == 0 ? 3 : v)
                                    / El::BigFloat(7 + seed + 2 * i));
      }
    return p;
  }

  // Newton divided differences: an independent evaluation of the
  // interpolant of the values f at the points x, at y.
  Boost_Float newton_interpolant(const std::vector<Boost_Float> &x,
                                 std::vector<Boost_Float> f,
                                 const Boost_Float &y)
  {
    const size_t K = x.size();
    for(size_t j = 1; j < K; ++j)
      for(size_t i = K - 1; i >= j; --i)
        f[i] = (f[i] - f[i - 1]) / (x[i] - x[i - j]);
    Boost_Float r = f[K - 1];
    for(size_t i = K - 1; i-- > 0;)
      r = r * (y - x[i]) + f[i];
    return r;
  }

  Boost_Float eval(const Polynomial &p, const Boost_Float &x)
  {
    Boost_Float r = 0;
    for(size_t i = p.coefficients.size(); i-- > 0;)
      r = r * x + to_Boost_Float(p.coefficients[i]);
    return r;
  }

  Polynomial_Vector_Matrix
  make_pvm(const Simple_Matrix<Polynomial_Vector> &polys,
           const Damped_Rational &prefactor,
           const std::optional<Damped_Rational> &reduced = std::nullopt)
  {
    return Polynomial_Vector_Matrix(polys, prefactor, reduced, std::nullopt,
                                    std::nullopt, std::nullopt, std::nullopt,
                                    std::nullopt);
  }

  // log2 of the largest |a - b| over entries (rows are orthonormal, so
  // entries are O(1))
  double log2_max_diff(const El::Matrix<El::BigFloat> &a,
                       const El::Matrix<El::BigFloat> &b)
  {
    REQUIRE(a.Height() == b.Height());
    REQUIRE(a.Width() == b.Width());
    El::BigFloat m = 0;
    for(int i = 0; i < a.Height(); ++i)
      for(int j = 0; j < a.Width(); ++j)
        m = std::max(m, El::Abs(a.CRef(i, j) - b.CRef(i, j)));
    return m == El::BigFloat(0)
             ? -1e9
             : static_cast<double>(log2(to_Boost_Float(m)));
  }

  // log2 of max |V V^T - 1|
  double log2_orthonormality_error(const El::Matrix<El::BigFloat> &v)
  {
    El::BigFloat m = 0;
    for(int i = 0; i < v.Height(); ++i)
      for(int j = 0; j < v.Height(); ++j)
        {
          El::BigFloat s = i == j ? -1 : 0;
          for(int k = 0; k < v.Width(); ++k)
            s += v.CRef(i, k) * v.CRef(j, k);
          m = std::max(m, El::Abs(s));
        }
    return m == El::BigFloat(0)
             ? -1e9
             : static_cast<double>(log2(to_Boost_Float(m)));
  }
}

TEST_CASE("two_set_encoding")
{
  if(El::mpi::Rank() != 0)
    return;
  const double precision = El::gmp::Precision();
  const Boost_Float pi = boost::math::constants::pi<Boost_Float>();

  SECTION("bases by recurrence = Cholesky bases on the default points")
  {
    // exp(-x) with two poles, degree 12: the moments' Cholesky is accurate
    // here, so both constructions must give the same sampled bases.
    const Damped_Rational prefactor{
      1, Boost_Float(1) / boost::math::constants::e<Boost_Float>(), {-1, -2}};
    const size_t degree = 12;
    const auto boost_points = sample_points(degree + 1, prefactor);
    const auto points = to_BigFloat_Vector(boost_points);
    const auto scalings
      = to_BigFloat_Vector(sample_scalings(boost_points, prefactor));
    const auto basis = bilinear_basis(points, scalings);
    const auto cholesky = sample_bilinear_basis(basis, points, scalings);
    const auto recurrence = sample_bilinear_basis_by_recurrence(
      {basis[0].size(), basis[1].size()}, points, scalings);
    for(size_t parity : {0, 1})
      {
        CAPTURE(parity);
        REQUIRE(log2_max_diff(cholesky[parity], recurrence[parity])
                < -(precision - 64));
        REQUIRE(log2_orthonormality_error(recurrence[parity])
                < -(precision - 32));
      }
    // Control: moving one point changes the recurrence's bases
    auto moved = points;
    moved[3] *= El::BigFloat("1.000001");
    const auto recurrence_moved = sample_bilinear_basis_by_recurrence(
      {basis[0].size(), basis[1].size()}, moved, scalings);
    REQUIRE(log2_max_diff(cholesky[0], recurrence_moved[0]) > -40);
  }

  SECTION("functional-matched without poles = Chebyshev points of x/(x+a)")
  {
    // g = (x + a)^(d+1): N_g = (2(d+1)/pi) atan sqrt(x/a), so the points
    // are a tan^2((2k-1) pi / (4(d+1))) and P_g = (d+1) sqrt(a) (x+a)^d / pi
    const Damped_Rational none{1, Boost_Float("0.2"), {}};
    const size_t num_points = 31;
    const El::BigFloat a = encoding_points_scale(none, El::BigFloat(0.5));
    REQUIRE(a == El::BigFloat(5));
    std::vector<El::BigFloat> fm_points, fm_scalings, ch_points, ch_scalings;
    functional_matched_encoding_points(none, num_points, a, fm_points,
                                       fm_scalings);
    chebyshev_encoding_points(none, num_points, a, ch_points, ch_scalings);
    const Boost_Float a_bf = to_Boost_Float(a);
    for(size_t k = 0; k < num_points; ++k)
      {
        CAPTURE(k);
        const Boost_Float x = to_Boost_Float(ch_points[k]);
        REQUIRE(rel_err(to_Boost_Float(fm_points[k]), x)
                < pow(Boost_Float(2), -(precision - 16)));
        const Boost_Float P_g
          = num_points * sqrt(a_bf) * pow(x + a_bf, int(num_points - 1)) / pi;
        REQUIRE(rel_err(to_Boost_Float(fm_scalings[k]), 1 / P_g)
                < pow(Boost_Float(2), -(precision - 16)));
      }
  }

  SECTION("functional-matched with poles: quadrature exact to degree d")
  {
    // d + 1 = 6 points, 5 poles, so g = prod (x + q_i) (x + a): with
    // c = {q_i, a} distinct, int_0^inf x^n dx / (sqrt(x) g(x))
    //   = pi sum_j (-c_j)^n / (sqrt(c_j) prod_{l != j} (c_l - c_j))
    // for n <= d (residues of z^(2n) / g(z^2) in the upper half-plane).
    const Damped_Rational poles{
      1, Boost_Float("0.2"), {Boost_Float("-0.3"), Boost_Float("-1.1"), Boost_Float("-2.7"),
                            Boost_Float(-6), Boost_Float(-15)}};
    const size_t num_points = 6;
    const El::BigFloat a = encoding_points_scale(poles, El::BigFloat(0.5));
    REQUIRE(to_Boost_Float(a) == Boost_Float("2.7") / 2);
    std::vector<El::BigFloat> points, scalings;
    functional_matched_encoding_points(poles, num_points, a, points,
                                       scalings);
    std::vector<Boost_Float> c{Boost_Float("0.3"), Boost_Float("1.1"),
                               Boost_Float("2.7"), Boost_Float(6),
                               Boost_Float(15), to_Boost_Float(a)};
    for(size_t n = 0; n < num_points; ++n)
      {
        Boost_Float exact = 0;
        for(size_t j = 0; j < c.size(); ++j)
          {
            Boost_Float den = sqrt(c[j]);
            for(size_t l = 0; l < c.size(); ++l)
              if(l != j)
                den *= c[l] - c[j];
            exact += pow(-c[j], int(n)) / den;
          }
        exact *= pi;
        Boost_Float sum = 0;
        for(size_t k = 0; k < num_points; ++k)
          sum += to_Boost_Float(scalings[k])
                 * pow(to_Boost_Float(points[k]), int(n));
        CAPTURE(n);
        REQUIRE(rel_err(sum, exact) < pow(Boost_Float(2), -(precision - 24)));
      }
    // Control: a scaling perturbed by 1e-6 breaks the exactness
    {
      auto perturbed = scalings;
      perturbed[2] *= El::BigFloat("1.000001");
      const int n = num_points - 1;
      Boost_Float exact = 0, sum = 0;
      for(size_t j = 0; j < c.size(); ++j)
        {
          Boost_Float den = sqrt(c[j]);
          for(size_t l = 0; l < c.size(); ++l)
            if(l != j)
              den *= c[l] - c[j];
          exact += pow(-c[j], n) / den;
        }
      for(size_t k = 0; k < num_points; ++k)
        sum += to_Boost_Float(perturbed[k])
               * pow(to_Boost_Float(points[k]), n);
      REQUIRE(rel_err(sum, pi * exact) > Boost_Float("1e-12"));
    }
  }

  SECTION("second step: rows = w_j times the interpolant at y_j")
  {
    // A 2x2 block in the reduced format: prefactor poles {-1,-2,-3,-20,-30},
    // reduced {-1,-2,-3}; degree 9, so 8 approximation points and the
    // program is the degree-7 interpolant of P / ((x+20)(x+30)).
    const Damped_Rational prefactor{
      1, Boost_Float("0.17"), {-1, -2, -3, -20, -30}};
    const Damped_Rational reduced{1, Boost_Float("0.17"), {-1, -2, -3}};
    const size_t degree = 9, vector_dim = 3;
    Simple_Matrix<Polynomial_Vector> polys(2, 2);
    for(int r = 0; r < 2; ++r)
      for(int s = 0; s < 2; ++s)
        {
          polys(r, s).clear();
          for(size_t n = 0; n < vector_dim; ++n)
            polys(r, s).push_back(
              test_polynomial(degree, 3 * std::min(r, s) + std::max(r, s)
                                        + 5 * static_cast<int>(n)));
        }
    auto m = make_pvm(polys, prefactor, reduced);
    const size_t K = m.sample_points.size();
    REQUIRE(K == 8);

    const Dual_Constraint_Group group_default(0, m);
    std::vector<El::BigFloat> points, scalings;
    functional_matched_encoding_points(
      m.reduced_prefactor, K,
      encoding_points_scale(m.reduced_prefactor, El::BigFloat(1)), points,
      scalings);
    // the points reach far beyond the approximation points
    REQUIRE(points.back() > 10 * m.sample_points.back());
    set_encoding_points(m, points, scalings);
    const Dual_Constraint_Group group(0, m);
    REQUIRE(group.num_points == K);
    REQUIRE(group.constraint_matrix.Height() == static_cast<El::Int>(3 * K));
    REQUIRE(group.constraint_matrix.Width() == vector_dim - 1);
    CAPTURE(group.encoding_digits_lost);
    REQUIRE(group.encoding_digits_lost > 0);
    REQUIRE(group.encoding_digits_lost < precision * std::log10(2.0) / 2);

    // Reference: Newton form of the interpolant of R at the approximation
    // points, from the polynomials' coefficients.
    const auto ref_row = [&](const int r, const int s, const size_t n,
                             const Boost_Float &y) {
      std::vector<Boost_Float> x, f;
      for(size_t k = 0; k < K; ++k)
        {
          const Boost_Float xk = to_Boost_Float(m.sample_points[k]);
          x.push_back(xk);
          f.push_back(eval(polys(r, s)[n], xk) / ((xk + 20) * (xk + 30)));
        }
      return newton_interpolant(x, f, y);
    };
    const Boost_Float tolerance
      = pow(Boost_Float(10), group.encoding_digits_lost)
        * pow(Boost_Float(2), -(precision - 24));
    const Boost_Float tolerance_default
      = pow(Boost_Float(2), -(precision - 24));
    size_t p = 0;
    for(int s = 0; s < 2; ++s)
      for(int r = 0; r <= s; ++r)
        for(size_t j = 0; j < K; ++j, ++p)
          {
            CAPTURE(r, s, j);
            const Boost_Float y = to_Boost_Float(points[j]);
            const Boost_Float w = to_Boost_Float(scalings[j]);
            Boost_Float row_max = 0;
            std::vector<Boost_Float> ref(vector_dim);
            for(size_t n = 0; n < vector_dim; ++n)
              {
                ref[n] = w * ref_row(r, s, n, y) * (n == 0 ? 1 : -1);
                row_max = std::max(row_max, Boost_Float(abs(ref[n])));
              }
            REQUIRE(abs(to_Boost_Float(group.constraint_constants[p]) - ref[0])
                    < tolerance * row_max);
            for(size_t n = 1; n < vector_dim; ++n)
              REQUIRE(abs(to_Boost_Float(group.constraint_matrix.CRef(p, n - 1))
                          - ref[n])
                      < tolerance * row_max);
            // The default encoding's rows are r_k R(x_k)
            const Boost_Float x = to_Boost_Float(m.sample_points[j]);
            const Boost_Float rs = to_Boost_Float(m.reduced_sample_scalings[j]);
            const Boost_Float ref_default
              = rs * eval(polys(r, s)[0], x) / ((x + 20) * (x + 30));
            REQUIRE(rel_err(to_Boost_Float(group_default.constraint_constants[p]),
                            ref_default)
                    < tolerance_default);
          }
    // Bases: recurrence at the encoding points, orthonormal rows
    const auto bases = sample_bilinear_basis_by_recurrence(
      {m.bilinear_basis[0].size(), m.bilinear_basis[1].size()}, points,
      scalings);
    for(size_t parity : {0, 1})
      {
        CAPTURE(parity);
        REQUIRE(group.bilinear_bases[parity].Height()
                == bases[parity].Height());
        REQUIRE(group.bilinear_bases[parity].Width() == static_cast<El::Int>(K));
        for(int i = 0; i < bases[parity].Height(); ++i)
          for(size_t k = 0; k < K; ++k)
            REQUIRE(group.bilinear_bases[parity].CRef(i, k)
                    == bases[parity].CRef(i, k));
        REQUIRE(log2_orthonormality_error(group.bilinear_bases[parity])
                < -(precision - 32));
      }

    // Control: a perturbed encoding point gives rows that fail the check
    {
      auto perturbed = points;
      perturbed[K - 1] *= El::BigFloat("1.000001");
      auto m2 = m;
      set_encoding_points(m2, perturbed, scalings);
      const Dual_Constraint_Group group2(0, m2);
      const Boost_Float y = to_Boost_Float(points[K - 1]);
      const Boost_Float w = to_Boost_Float(scalings[K - 1]);
      const Boost_Float ref = w * ref_row(0, 0, 0, y);
      REQUIRE(abs(to_Boost_Float(group2.constraint_constants[K - 1]) - ref)
              > 1000 * tolerance * abs(ref));
    }
  }

  SECTION("without a reduced prefactor the two-set rows are w_j P(y_j)")
  {
    // The interpolant of a degree-d polynomial at d + 1 points is itself.
    const Damped_Rational prefactor{1, Boost_Float("0.3"), {-0.5, -4}};
    const size_t degree = 14;
    Simple_Matrix<Polynomial_Vector> polys(1, 1);
    polys(0, 0) = {test_polynomial(degree, 1), test_polynomial(degree, 2)};
    auto m = make_pvm(polys, prefactor);
    const size_t K = m.sample_points.size();
    REQUIRE(K == degree + 1);
    std::vector<El::BigFloat> points, scalings;
    chebyshev_encoding_points(
      m.reduced_prefactor, K,
      encoding_points_scale(m.reduced_prefactor, El::BigFloat(0.5)), points,
      scalings);
    set_encoding_points(m, points, scalings);
    const Dual_Constraint_Group group(0, m);
    const Boost_Float tolerance
      = pow(Boost_Float(10), group.encoding_digits_lost)
        * pow(Boost_Float(2), -(precision - 24));
    for(size_t j = 0; j < K; ++j)
      {
        CAPTURE(j);
        const Boost_Float y = to_Boost_Float(points[j]);
        const Boost_Float w = to_Boost_Float(scalings[j]);
        const Boost_Float c = w * eval(polys(0, 0)[0], y);
        const Boost_Float b = -w * eval(polys(0, 0)[1], y);
        const Boost_Float row_max = std::max(Boost_Float(abs(c)), Boost_Float(abs(b)));
        REQUIRE(abs(to_Boost_Float(group.constraint_constants[j]) - c)
                < tolerance * row_max);
        REQUIRE(abs(to_Boost_Float(group.constraint_matrix.CRef(j, 0)) - b)
                < tolerance * row_max);
      }
  }

  SECTION("degree 1 (two points)")
  {
    // The smallest block the two-set step re-encodes.  Without a reduced
    // prefactor (P of degree 1) the rows are w_j P(y_j); with one (P of
    // degree 2, one pole removed), w_j times the line through
    // R = P / (x + 4) at the two approximation points.
    const Damped_Rational prefactor{1, Boost_Float("0.3"), {-0.5, -4}};
    const Damped_Rational reduced{1, Boost_Float("0.3"), {-0.5}};
    for(const bool with_reduced : {false, true})
      for(const bool chebyshev : {false, true})
        {
          CAPTURE(with_reduced, chebyshev);
          Simple_Matrix<Polynomial_Vector> polys(1, 1);
          // with the reduced prefactor one pole is removed, so P has degree
          // 2 and R degree 1
          const size_t degree = with_reduced ? 2 : 1;
          polys(0, 0)
            = {test_polynomial(degree, 1), test_polynomial(degree, 2)};
          auto m = with_reduced ? make_pvm(polys, prefactor, reduced)
                                : make_pvm(polys, prefactor);
          const size_t K = m.sample_points.size();
          REQUIRE(K == 2);
          const auto a
            = encoding_points_scale(m.reduced_prefactor, El::BigFloat(1));
          std::vector<El::BigFloat> points, scalings;
          if(chebyshev)
            chebyshev_encoding_points(m.reduced_prefactor, K, a, points,
                                      scalings);
          else
            functional_matched_encoding_points(m.reduced_prefactor, K, a,
                                               points, scalings);
          REQUIRE(points.size() == K);
          set_encoding_points(m, points, scalings);
          const Dual_Constraint_Group group(0, m);
          REQUIRE(group.num_points == K);
          const Boost_Float tolerance
            = pow(Boost_Float(10), group.encoding_digits_lost)
              * pow(Boost_Float(2), -(precision - 24));
          // R at the approximation points, and its line at y
          const auto R = [&](const size_t n, const Boost_Float &x) {
            Boost_Float r = eval(polys(0, 0)[n], x);
            if(with_reduced)
              r /= x + 4;
            return r;
          };
          const Boost_Float x0 = to_Boost_Float(m.sample_points[0]),
                            x1 = to_Boost_Float(m.sample_points[1]);
          for(size_t j = 0; j < K; ++j)
            {
              CAPTURE(j);
              const Boost_Float y = to_Boost_Float(points[j]);
              const Boost_Float w = to_Boost_Float(scalings[j]);
              const auto line = [&](const size_t n) {
                return R(n, x0) + (R(n, x1) - R(n, x0)) * (y - x0) / (x1 - x0);
              };
              const Boost_Float c = w * line(0), b = -w * line(1);
              const Boost_Float row_max
                = std::max(Boost_Float(abs(c)), Boost_Float(abs(b)));
              REQUIRE(abs(to_Boost_Float(group.constraint_constants[j]) - c)
                      < tolerance * row_max);
              REQUIRE(abs(to_Boost_Float(group.constraint_matrix.CRef(j, 0))
                          - b)
                      < tolerance * row_max);
            }
          for(size_t parity : {0, 1})
            {
              CAPTURE(parity);
              REQUIRE(group.bilinear_bases[parity].Width()
                      == static_cast<El::Int>(K));
              REQUIRE(group.bilinear_bases[parity].Height() == 1);
              REQUIRE(log2_orthonormality_error(group.bilinear_bases[parity])
                      < -(precision - 32));
            }
        }
  }

  SECTION("degree 0 (one point)")
  {
    // A constant block: one point, an even basis of one element and an empty
    // odd basis.  (pmp2sdp leaves such blocks on the default path; the
    // second step must still handle them.)  Moving its point to y with
    // scaling w gives the row w R(x_0): R is constant.
    const Damped_Rational prefactor{1, Boost_Float("0.3"), {}};
    Simple_Matrix<Polynomial_Vector> polys(1, 1);
    polys(0, 0) = {test_polynomial(0, 1), test_polynomial(0, 2)};
    auto m = make_pvm(polys, prefactor);
    REQUIRE(m.sample_points.size() == 1);
    const Dual_Constraint_Group group_default(0, m);
    const El::BigFloat y("2.5"), w("0.75");
    set_encoding_points(m, {y}, {w});
    const Dual_Constraint_Group group(0, m);
    REQUIRE(group.num_points == 1);
    REQUIRE(group.encoding_digits_lost < 1e-10);
    const El::BigFloat sigma = m.reduced_sample_scalings[0];
    REQUIRE(rel_err(to_Boost_Float(group.constraint_constants[0]),
                    to_Boost_Float(w * group_default.constraint_constants[0]
                                   / sigma))
            < pow(Boost_Float(2), -(precision - 8)));
    REQUIRE(rel_err(to_Boost_Float(group.constraint_matrix.CRef(0, 0)),
                    to_Boost_Float(w * group_default.constraint_matrix.CRef(0, 0)
                                   / sigma))
            < pow(Boost_Float(2), -(precision - 8)));
    REQUIRE(group.bilinear_bases[0].Height() == 1);
    REQUIRE(group.bilinear_bases[0].Width() == 1);
    REQUIRE(rel_err(to_Boost_Float(group.bilinear_bases[0].CRef(0, 0)), 1)
            < pow(Boost_Float(2), -(precision - 8)));
    REQUIRE(group.bilinear_bases[1].Height() == 0);
    REQUIRE(group.bilinear_bases[1].Height()
            == group_default.bilinear_bases[1].Height());
  }

  SECTION("encoding points = approximation points reproduce the default")
  {
    const Damped_Rational prefactor{1, Boost_Float("0.3"), {-0.5, -4, -9}};
    const Damped_Rational reduced{1, Boost_Float("0.3"), {-0.5}};
    Simple_Matrix<Polynomial_Vector> polys(1, 1);
    polys(0, 0) = {test_polynomial(12, 3), test_polynomial(12, 4)};
    auto m = make_pvm(polys, prefactor, reduced);
    const Dual_Constraint_Group group_default(0, m);
    set_encoding_points(m, m.sample_points, m.reduced_sample_scalings);
    const Dual_Constraint_Group group(0, m);
    REQUIRE(group.encoding_digits_lost == 0);
    for(size_t j = 0; j < m.sample_points.size(); ++j)
      {
        CAPTURE(j);
        REQUIRE(rel_err(to_Boost_Float(group.constraint_constants[j]),
                        to_Boost_Float(group_default.constraint_constants[j]))
                < pow(Boost_Float(2), -(precision - 8)));
        REQUIRE(rel_err(to_Boost_Float(group.constraint_matrix.CRef(j, 0)),
                        to_Boost_Float(group_default.constraint_matrix.CRef(j, 0)))
                < pow(Boost_Float(2), -(precision - 8)));
      }
  }

  SECTION("binary output at a lower precision: El's limb count")
  {
    // pmp2sdp --encodingPoints writes BigFloats with the output precision's
    // limb count; the formula must be El's own.
    REQUIRE(sdpb_serialization::num_limbs_for_precision(El::gmp::Precision())
            == static_cast<int>(El::BigFloat().NumLimbs()));
    REQUIRE(sdpb_serialization::output_num_limbs() == 0);
  }

  SECTION("set_encoding_points rejects bad input")
  {
    const Damped_Rational prefactor{1, Boost_Float("0.3"), {-0.5}};
    Simple_Matrix<Polynomial_Vector> polys(1, 1);
    polys(0, 0) = {test_polynomial(4, 1), test_polynomial(4, 2)};
    auto m = make_pvm(polys, prefactor);
    std::vector<El::BigFloat> pts{1, 2, 3, 4, 5}, sc{1, 1, 1, 1, 1};
    REQUIRE_THROWS(set_encoding_points(m, {1, 2, 3, 4}, {1, 1, 1, 1}));
    REQUIRE_THROWS(set_encoding_points(m, {1, 2, 2, 4, 5}, sc));
    REQUIRE_THROWS(set_encoding_points(m, pts, {1, 1, 0, 1, 1}));
    REQUIRE_THROWS(set_encoding_points(m, {-1, 2, 3, 4, 5}, sc));
    REQUIRE_NOTHROW(set_encoding_points(m, pts, sc));
  }
}
