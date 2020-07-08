#include <El.hpp>

template<typename Real>
void RandomFeasibleLP( El::Int m, El::Int n, El::Int k )
{
  El::Output("Testing with ",El::TypeName<Real>());
  // Create random (primal feasible) inputs for the primal/dual problem
  //    arginf_{x,s} { c^T x | A x = b, G x + s = h, s >= 0 }
  //    argsup_{y,z} { -b^T y - h^T z | A^T y + G^T z + c = 0, z >= 0 }.

  // xFeas and sFeas are only used for problem generation
  El::Matrix<Real> xFeas, sFeas;
  El::Uniform( xFeas, n, 1 ); // Sample over B_1(0)
  El::Uniform( sFeas, k, 1, Real(1), Real(1) ); // Sample over B_1(1)

  El::Matrix<Real> A, G, b, c, h;
  El::Uniform( A, m, n );
  El::Uniform( G, k, n );
  El::Gemv( El::NORMAL, Real(1), A, xFeas, b );
  El::Gemv( El::NORMAL, Real(1), G, xFeas, h );
  h += sFeas;
  El::Uniform( c, n, 1 );

  // Solve the primal/dual Linear Program with the default options
  El::Matrix<Real> x, y, z, s;
  El::Timer timer;
  timer.Start();
  El::lp::affine::Ctrl<Real> ctrl;
  ctrl.mehrotraCtrl.print=true;
  // ctrl.mehrotraCtrl.maxIts=10;
  El::LP( A, G, b, c, h, x, y, z, s, ctrl );
  // El::LP( A, G, b, c, h, x, y, z, s );
  El::Output("Primal-dual LP took ",timer.Stop()," seconds");

  // Print the primal and dual objective values
  const Real primal = El::Dot(c,x);
  const Real dual = -El::Dot(b,y) - El::Dot(h,z);
  const Real relGap = El::Abs(primal-dual) / El::Max(El::Abs(dual),Real(1));
  El::Output("c^T x = ",primal);
  El::Output("-b^T y - h^T z = ",dual);
  El::Output("|gap| / max( |dual|, 1 ) = ",relGap);

  // Print the relative primal feasibility residual,
  //   || A x - b ||_2 / max( || b ||_2, 1 ).
  El::Matrix<Real> rPrimal;
  El::Gemv( El::NORMAL, Real(1), A, x, rPrimal );
  rPrimal -= b;
  const Real bFrob = El::FrobeniusNorm( b );
  const Real rPrimalFrob = El::FrobeniusNorm( rPrimal );
  const Real primalRelResid = rPrimalFrob / El::Max( bFrob, Real(1) );
  El::Output("|| A x - b ||_2 / || b ||_2 = ",primalRelResid);

  // Print the relative dual feasiability residual,
  //   || A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ).
  El::Matrix<Real> rDual;
  El::Gemv( El::TRANSPOSE, Real(1), A, y, rDual );
  El::Gemv( El::TRANSPOSE, Real(1), G, z, Real(1), rDual );
  rDual += c;
  const Real cFrob = El::FrobeniusNorm( c );
  const Real rDualFrob = El::FrobeniusNorm( rDual );
  const Real dualRelResid = rDualFrob / El::Max( cFrob, Real(1) );
  El::Output
  ("|| A^T y + G^T z + c ||_2 / max( || c ||_2, 1 ) = ",dualRelResid);
  El::Output("");
}

int main(int argc, char **argv)
{
  El::Environment env(argc, argv);

  // const El::Int m = El::Input("--m","height of A",70);
  // const El::Int n = El::Input("--n","width of A",80);
  // const El::Int k = El::Input("--k","height of G",90);
  // El::ProcessInput();

  // RandomFeasibleLP<float>( m, n, k );
  // RandomFeasibleLP<double>( m, n, k );
  El::gmp::SetPrecision(256);
  // RandomFeasibleLP<El::BigFloat>( m, n, k );


  
  {
    const int64_t n(3), m(1);
    El::Matrix<El::BigFloat> b(m,1);
    El::Matrix<El::BigFloat> A(m,n), c(n,1);

    c(0)=1;
    c(1)=-1;
    c(2)=0;

    A(0,0)="-8.34333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e6";
    A(0,1)=-A(0,0);
    A(0,2)=1;

    b(0)="1.00000001e8";
  
    // Solve the primal/dual Linear Program with the default options
    El::Matrix<El::BigFloat> x, y, z;
    El::lp::direct::Ctrl<El::BigFloat> ctrl(false);
    // ctrl.mehrotraCtrl.print=true;
    // ctrl.mehrotraCtrl.maxIts=10;
    El::LP( A, b, c, x, y, z, ctrl );
    std::cout.precision(1024/4);
    std::cout << (x(0)-x(1)) << "\n"
              // << x(0) << "\n"
              // << x(1) << "\n"
              // << x(2) << "\n"
      ;
    // El::Print(A,"A");
    // std::cout << "\n";
  }

  {
    const int64_t n(4), m(2);
    El::Matrix<El::BigFloat> b(m,1);
    El::Matrix<El::BigFloat> A(m,n), c(n,1);

    c(0)=1;
    c(1)=-1;
    c(2)=0;
    c(3)=0;

    A(0,0)="-8.34333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e6";
    A(0,1)=-A(0,0);
    A(0,2)=1;
    A(0,3)=0;

    b(0)="1.00000001e8";

    A(1,0)="-2.0883681393738418786521971610157869999564952571145363799181069183090349103427762822892226369869922745782002357632815274877025509325687194149604661156665774994079238964936547103062281911314159394645950806565436291836646e6";
    A(1,1)=-A(1,0);
    A(1,2)=0;
    A(1,3)=1;

    b(1)="2.500041817188193259583610690651477602462222452320541114185075697450625764505571202951722797034907521311172640839539068960526921868899591722427115600132424981324576028654978230857101759743577309108511615671925085573048409e7";
    
    // Solve the primal/dual Linear Program with the default options
    El::Matrix<El::BigFloat> x, y, z;
    El::lp::direct::Ctrl<El::BigFloat> ctrl(false);
    // ctrl.mehrotraCtrl.print=true;
    // ctrl.mehrotraCtrl.maxIts=10;
    El::LP( A, b, c, x, y, z, ctrl );
    std::cout.precision(1024/4);
    std::cout << (x(0)-x(1)) << "\n"
              // << x(0) << "\n"
              // << x(1) << "\n"
              // << x(2) << "\n"
              // << x(3) << "\n"
      ;
    // El::Print(A,"A");
    // std::cout << "\n";
    // El::Print(b,"b");
    // std::cout << "\n";
  }

  {
    const int64_t n(5), m(3);
    El::Matrix<El::BigFloat> b(m,1);
    El::Matrix<El::BigFloat> A(m,n), c(n,1);

    c(0)=1;
    c(1)=-1;
    c(2)=0;
    c(3)=0;
    c(4)=0;

    A(0,0)="-8.34333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333e6";
    A(0,1)=-A(0,0);
    A(0,2)=1;
    A(0,3)=0;
    A(0,4)=0;

    b(0)="1.00000001e8";

    A(1,0)="-2.0883681393738418786521971610157869999564952571145363799181069183090349103427762822892226369869922745782002357632815274877025509325687194149604661156665774994079238964936547103062281911314159394645950806565436291836646e6";
    A(1,1)=-A(1,0);
    A(1,2)=0;
    A(1,3)=1;
    A(1,4)=0;

    b(1)="2.500041817188193259583610690651477602462222452320541114185075697450625764505571202951722797034907521311172640839539068960526921868899591722427115600132424981324576028654978230857101759743577309108511615671925085573048409e7";

    A(2,0)="-523359.49039448594574060530655172286187090378437496537397484654111114679353770021695730603264685384461177908057073681681768485010274393686753539800942415710439431439847258920370816312326141491421457729752660492298893732";
    A(2,1)=-A(2,0);
    A(2,2)=0;
    A(2,3)=0;
    A(2,4)=1;

    b(2)="6.2503141332235077264522964494075174019801929035358155014826774648939365836286594330371275920790232914052411157897491638288477145045802473293978949526674980879470479960005412142141993037669145499508595145049713993404956e6";


    
    // Solve the primal/dual Linear Program with the default options
    El::Matrix<El::BigFloat> x, y, z;
    El::lp::direct::Ctrl<El::BigFloat> ctrl(false);
    // ctrl.mehrotraCtrl.print=true;
    // ctrl.mehrotraCtrl.maxIts=10;
    El::LP( A, b, c, x, y, z, ctrl );
    std::cout.precision(1024/4);
    std::cout << (x(0)-x(1)) << "\n"
              // << x(0) << "\n"
              // << x(1) << "\n"
              // << x(2) << "\n"
              // << x(3) << "\n"
      ;
    // El::Print(A,"A");
    // std::cout << "\n";
    // El::Print(b,"b");
    // std::cout << "\n";
  }
  
}
