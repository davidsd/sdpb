#include "pmp/PMP_Info.hxx"
#include "spectrum/Zero_Test.hxx"
#include "spectrum/compute_spectrum/zero_depth.hxx"
#include "sdpb_util/Timers/Timers.hxx"

#include <catch2/catch_amalgamated.hpp>

#include <El.hpp>
#include <vector>

// Defined in src/spectrum/compute_spectrum/find_zeros.cxx
std::vector<El::BigFloat>
find_zeros(const El::Matrix<El::BigFloat> &c_minus_By_block,
           const PVM_Info &pvm, const Zero_Test &zero_test,
           const Boost_Float &threshold, const El::BigFloat &max_zero,
           const El::BigFloat &min_zero_distance, Timers &timers);

// Tests for the zero tests of spectrum.
// Toy functionals have exactly known zeros:
//   P(x) = C prod_k (x + k/2) (x + delta) prod_i [(x - a_i)^2 + w_i^2] prod_j [(x - b_j)^2 + s_j^2]
// True zeros are a_i with tiny widths w_i (and x=0 if delta is tiny);
// b_j with sizeable s_j are minima that are not zeros.

namespace
{
  Boost_Polynomial quadratic_factor(const Boost_Float &a, const Boost_Float &w)
  {
    return Boost_Polynomial({a * a + w * w, -2 * a, Boost_Float(1)});
  }

  struct Toy
  {
    // Polynomial part (diagonal entries of the matrix)
    std::vector<Boost_Polynomial> diagonal;
    Damped_Rational prefactor;
    // Matrix = R^T diag(...) R, with unimodular R, so det = product of diagonal
    bool mix = true;
  };

  Damped_Rational exp_decay(const Boost_Float &base, const size_t num_poles)
  {
    Damped_Rational result;
    result.constant = 1;
    result.base = base;
    for(size_t k = 0; k < num_poles; ++k)
      result.poles.push_back(-Boost_Float(k) - Boost_Float(1) / 2);
    return result;
  }

  // Build PVM_Info and (c-B.y) as SDPB writes them:
  // c_minus_By[r,s,k] = M_{r,s}(x_k) * chi(x_k) * scale_k,
  // reduced_sample_scalings[k] = chi(x_k) * scale_k.
  std::pair<PVM_Info, El::Matrix<El::BigFloat>>
  make_block(const Toy &toy, const Boost_Float &x_max)
  {
    const size_t dim = toy.diagonal.size();
    size_t degree = 0;
    for(auto &p : toy.diagonal)
      degree = std::max(degree, p.degree());
    const size_t num_points = degree + 1;

    PVM_Info pvm;
    pvm.block_index = 0;
    pvm.block_path = "toy";
    pvm.dim = dim;
    pvm.prefactor = toy.prefactor;
    pvm.reduced_prefactor = toy.prefactor;
    for(size_t k = 0; k < num_points; ++k)
      {
        // points accumulate near x=0, as Laguerre-type sample points do
        const Boost_Float t = Boost_Float(k) / (num_points - 1);
        const Boost_Float x = x_max * t * t + x_max * t / 1000;
        pvm.sample_points.push_back(to_BigFloat(x));
        const Boost_Float scaling = toy.prefactor.evaluate(x);
        pvm.sample_scalings.push_back(to_BigFloat(scaling));
        pvm.reduced_sample_scalings.push_back(to_BigFloat(scaling));
      }

    El::Matrix<El::BigFloat> c_minus_By(dim * (dim + 1) / 2 * num_points, 1);
    int index = 0;
    for(size_t i = 0; i < dim; ++i)
      for(size_t j = 0; j <= i; ++j)
        for(size_t k = 0; k < num_points; ++k)
          {
            const Boost_Float x = to_Boost_Float(pvm.sample_points.at(k));
            // M = R^T D R with R = 1 + (upper shift), i.e.
            // M_ij = sum_l R_li D_l R_lj, R_li = delta_li + delta_{l+1,i}
            Boost_Float value = 0;
            for(size_t l = 0; l < dim; ++l)
              {
                const auto R = [&](size_t row, size_t col) -> Boost_Float {
                  if(row == col)
                    return 1;
                  if(toy.mix && row + 1 == col)
                    return 1;
                  return 0;
                };
                value
                  += R(l, i) * toy.diagonal.at(l).evaluate(x) * R(l, j);
              }
            c_minus_By(index, 0) = to_BigFloat(
              value * to_Boost_Float(pvm.reduced_sample_scalings.at(k)));
            ++index;
          }
    return {pvm, c_minus_By};
  }

  Boost_Polynomial background(const size_t q, const Boost_Float &constant)
  {
    Boost_Polynomial result({constant});
    for(size_t k = 1; k <= q; ++k)
      result *= Boost_Polynomial({Boost_Float(k) / 2, Boost_Float(1)});
    return result;
  }

  std::vector<El::BigFloat>
  run(const Toy &toy, const Boost_Float &x_max, const Zero_Test &test,
      const Boost_Float &threshold = Boost_Float("1e-10"))
  {
    auto [pvm, c_minus_By] = make_block(toy, x_max);
    Timers timers;
    return find_zeros(c_minus_By, pvm, test, threshold, 0, 0, timers);
  }

  // Check that zeros are close to expected values (relative tolerance)
  void require_zeros(const std::vector<El::BigFloat> &zeros,
                     const std::vector<double> &expected)
  {
    CAPTURE(zeros);
    CAPTURE(expected);
    REQUIRE(zeros.size() == expected.size());
    for(size_t i = 0; i < zeros.size(); ++i)
      {
        CAPTURE(i);
        REQUIRE(El::Abs(zeros.at(i) - El::BigFloat(expected.at(i)))
                <= El::BigFloat(1e-8) * (1 + El::BigFloat(expected.at(i))));
      }
  }
}

TEST_CASE("spectrum_zero_depth")
{
  // 1x1 matrices: M = P, so P'/P and P''/P are exact
  const auto depth_1x1 = [](const Boost_Polynomial &P, const Boost_Float &x,
                            const Boost_Float &L, const bool boundary) {
    El::Matrix<El::BigFloat> M(1, 1), dM(1, 1), d2M(1, 1);
    M(0, 0) = to_BigFloat(P.evaluate(x));
    dM(0, 0) = to_BigFloat(P.prime().evaluate(x));
    d2M(0, 0) = to_BigFloat(P.prime().prime().evaluate(x));
    return zero_local_depth(M, dM, d2M, to_BigFloat(L), boundary);
  };
  const Boost_Float a = 3, q = Boost_Float("1e40");
  const Boost_Float w = Boost_Float("1e-20");

  SECTION("quadratic dip: depth = (w/L)^2")
  {
    const auto P = quadratic_factor(a, w) * Boost_Polynomial({q});
    const Boost_Float L = Boost_Float(1) / 2;
    const auto depth = depth_1x1(P, a, L, false);
    REQUIRE(abs(depth / (w * w / (L * L)) - 1) < Boost_Float("1e-30"));
  }

  SECTION("not a minimum: the first-order term rejects it")
  {
    // P = q (x + 1)^2 at x = 3: P'/P = 1/2, P''/P = 1/8; L = 1
    const auto P = quadratic_factor(-1, 0) * Boost_Polynomial({q});
    const auto depth = depth_1x1(P, a, 1, false);
    REQUIRE(abs(depth - 2) < Boost_Float("1e-30"));
  }

  SECTION("boundary: simple zero just left of x=0")
  {
    const Boost_Float delta = Boost_Float("1e-30"), L = 2;
    const Boost_Polynomial P({q * delta, q});
    REQUIRE(abs(depth_1x1(P, 0, L, true) / (delta / L) - 1)
            < Boost_Float("1e-30"));
    // decreasing at x=0: not a zero
    REQUIRE(isinf(depth_1x1(Boost_Polynomial({q, -q}), 0, L, true)));
  }

  SECTION("negative value is a zero")
  {
    const auto P = quadratic_factor(a, 0) + Boost_Polynomial({-w});
    REQUIRE(depth_1x1(P, a, 1, false) < 0);
  }

  SECTION("rank-2 dip: depth of the narrowest eigenvalue, not the product")
  {
    // M = diag((x-a)^2 + w1^2, (x-a)^2 + w2^2), mixed by a unimodular R
    const Boost_Float w1("1e-25"), w2("1e-12"), L = 1;
    El::Matrix<El::BigFloat> D(2, 2), dD(2, 2), d2D(2, 2), R(2, 2);
    El::Zero(D);
    El::Zero(dD);
    El::Zero(d2D);
    D(0, 0) = to_BigFloat(w1 * w1);
    D(1, 1) = to_BigFloat(w2 * w2);
    d2D(0, 0) = d2D(1, 1) = 2;
    R(0, 0) = R(0, 1) = R(1, 1) = 1;
    R(1, 0) = 0;
    const auto mix = [&](const El::Matrix<El::BigFloat> &X) {
      El::Matrix<El::BigFloat> T, result;
      El::Gemm(El::TRANSPOSE, El::NORMAL, El::BigFloat(1), R, X, T);
      El::Gemm(El::NORMAL, El::NORMAL, El::BigFloat(1), T, R, result);
      return result;
    };
    const auto depth = zero_local_depth(mix(D), mix(dD), mix(d2D), 1, false);
    // 1 / (1/w1^2 + 1/w2^2) ~ w1^2
    const Boost_Float expected = 1 / (1 / (w1 * w1) + 1 / (w2 * w2));
    REQUIRE(abs(depth / expected - 1) < Boost_Float("1e-30"));
  }

  SECTION("local length")
  {
    const std::vector<El::BigFloat> crit{El::BigFloat(1), El::BigFloat(2),
                                         El::BigFloat(5)};
    const El::BigFloat zero(0), fallback(100);
    // nearest critical point
    REQUIRE(zero_local_length(El::BigFloat(2), false, crit, zero, fallback)
            == El::BigFloat(1));
    // no critical point on the left: the boundary counts
    const std::vector<El::BigFloat> crit2{El::BigFloat("0.1"),
                                          El::BigFloat(5)};
    REQUIRE(zero_local_length(El::BigFloat("0.1"), false, crit2, zero,
                              fallback)
            == El::BigFloat("0.1"));
    // boundary candidate: first critical point
    REQUIRE(zero_local_length(zero, true, crit, zero, fallback)
            == El::BigFloat(1));
    // nothing else: fallback
    REQUIRE(zero_local_length(zero, true, {}, zero, fallback) == fallback);
    // monotone P: a spurious far root of P' must not set L
    const std::vector<El::BigFloat> crit3{El::BigFloat("1.6e44")};
    REQUIRE(zero_local_length(zero, true, crit3, zero, fallback) == fallback);
    REQUIRE(zero_local_length(El::BigFloat(3), false, crit3, zero, fallback)
            == El::BigFloat(3));
  }
}

TEST_CASE("spectrum_find_zeros")
{
  const Boost_Float r = 3 - 2 * sqrt(Boost_Float(2));

  SECTION("zeros and wiggles, no prefactor")
  {
    Toy toy;
    toy.prefactor = exp_decay(1, 0);
    auto P = background(0, 1);
    P *= quadratic_factor(1, Boost_Float("1e-25"));
    P *= quadratic_factor(4, Boost_Float("1e-25"));
    P *= quadratic_factor(Boost_Float("2.2"), Boost_Float("0.5"));
    toy.diagonal = {P};
    for(auto test : {Zero_Test::curvature, Zero_Test::geometric_mean})
      {
        CAPTURE(to_string(test));
        require_zeros(run(toy, 10, test), {1, 4});
      }
  }

  SECTION("wiggle in the exponential tail (false zero of geometricMean)")
  {
    // det(c-B.y) decays like r^(3x); the wiggle near x=118 is not a zero,
    // but it is ~1e-16 below the geometric mean of its distant neighbours.
    Toy toy;
    toy.prefactor = exp_decay(r, 20);
    auto P = background(20, 1);
    P *= quadratic_factor(1, Boost_Float("1e-25"));
    P *= quadratic_factor(50, Boost_Float("1e-20"));
    P *= quadratic_factor(225, Boost_Float("1e-18"));
    P *= quadratic_factor(120, 4);
    toy.diagonal = {P, background(0, 1), background(0, 1)};
    toy.mix = false;
    require_zeros(run(toy, 300, Zero_Test::curvature), {1, 50, 225});
    // Documents the failure of the old test
    REQUIRE(run(toy, 300, Zero_Test::geometric_mean).size() == 4);
  }

  SECTION("J-shaped zero")
  {
    Toy toy;
    toy.prefactor = exp_decay(r, 0);
    auto P = background(5, 1);
    P *= quadratic_factor(3, Boost_Float("1e-25"));
    P *= quadratic_factor(52, Boost_Float("1e-22"));
    P *= quadratic_factor(300, Boost_Float("1e-20"));
    toy.diagonal = {P, background(0, 1), background(0, 1)};
    toy.mix = false;
    require_zeros(run(toy, 400, Zero_Test::curvature), {3, 52, 300});
  }

  SECTION("boundary zero, rank-2 zero, 2x2 matrix")
  {
    Toy toy;
    toy.prefactor = exp_decay(r, 5);
    // diagonal entries: one with a boundary zero and a zero at 3,
    // the other with a zero at 3 (rank-2 zero) and at 10.
    auto P1 = Boost_Polynomial({Boost_Float("1e-25"), Boost_Float(1)});
    P1 *= quadratic_factor(3, Boost_Float("1e-22"));
    P1 *= quadratic_factor(Boost_Float("6.5"), 2);
    auto P2 = quadratic_factor(3, Boost_Float("1e-15"));
    P2 *= quadratic_factor(10, Boost_Float("1e-20"));
    toy.diagonal = {P1, P2};
    require_zeros(run(toy, 20, Zero_Test::curvature), {0, 3, 10});
  }

  SECTION("scale invariance: x -> 1000 x, P -> 1e100 P")
  {
    Toy toy;
    toy.prefactor = exp_decay(1, 0);
    auto P = background(0, Boost_Float("1e100"));
    P *= quadratic_factor(1000, Boost_Float("1e-22"));
    P *= quadratic_factor(4000, Boost_Float("1e-22"));
    P *= quadratic_factor(2200, 500);
    toy.diagonal = {P};
    require_zeros(run(toy, 10000, Zero_Test::curvature), {1000, 4000});
  }
}
