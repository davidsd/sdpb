#include <stdlib.h> 
#include "../../../dynamical_solve.hxx"
#include <iostream>
#include <El.hpp>
using namespace std;

void newton_search1d (const El::BigFloat &thresh, const int &it, const El::Matrix<El::Base<El::BigFloat>> &w, El::BigFloat &lambda, const El::Matrix<El::BigFloat> &e_dot_grad, const El::BigFloat &radius, const int &n_external_parameters){
   if (it == 0){
     throw std::invalid_argument( "Newton1dMaxIterationsExceeded");
   }
   else{
     El::BigFloat f(-radius * radius);
     El::BigFloat fprime(0) ; 
     for(int i = 0; i < n_external_parameters; i++){
       //cout << "for loop"<<endl;
       f += (e_dot_grad(i) *e_dot_grad(i))/((w(i) + lambda)*(w(i) + lambda)); 
       fprime -= 2* (e_dot_grad(i)*e_dot_grad(i))/((w(i) + lambda)*(w(i) + lambda)*(w(i) + lambda));
       //cout << "it: " << it << endl;
     } 
     El::BigFloat step(- f/fprime);
     //cout << "step: "<< step << endl;
     if (El::Abs(step) < thresh){
       //cout <<"smaller? " <<  "it: " << it << endl;
       return ;
     }
     else{
       lambda += step;
       newton_search1d (thresh, it-1, w, lambda, e_dot_grad, radius, n_external_parameters);
     }
   }
}


void trust_region_boundary_step (const El::BigFloat &thresh, const int &max_it, const El::BigFloat &radius, 
                                 El::Matrix<El::BigFloat> &A, const El::Matrix<El::BigFloat> &grad, El::Matrix<El::BigFloat> &step){ 
  //cout<<"start of function" << endl;
  int n_external_parameters = grad.Height();
  typedef El::Base<El::BigFloat> Real; 
  El::Matrix<Real> w(n_external_parameters,1);
  El::Zero(w);
  El::Matrix<El::BigFloat> Q(A);
  //cout << "start solving " << endl;
  //El::Matrix<El::BigFloat> householderScalars;
  //El::HermitianTridiag( El::UPPER, A, householderScalars );
  //auto d = El::GetRealPartOfDiagonal(A);
  //auto dSub = El::GetDiagonal( A, 1 );
  //El::HermitianTridiagEig( d, dSub, w, Q, ctrl.tridiagEigCtrl );
  //El::herm_tridiag::ApplyQ( El::LEFT, El::UPPER, El::NORMAL, A, householderScalars, Q );
  //El::apply_packed_reflectors::LLVBUnblocked(El::CONJUGATED, 1, Q, householderScalars, A);
  El::HermitianEig(El::LOWER, A, w, Q);//,ctrl);
  El::Print( w, "\neigenvalues:");
  El::Print( Q, "\neigenvectors:");
  //El::Print( Q(El::ALL, 0), "eigenvectors:" );
  //El::Print( Q(El::ALL, 1), "eigenvectors:" );
  El::BigFloat lambda0(0);   
  El::Matrix<El::BigFloat> e_dot_grad(grad);
  for (int i = 0; i < n_external_parameters; i++){
    e_dot_grad(i) = El::Abs(El::Dot(Q(El::ALL, i),grad));
    El::BigFloat temp (e_dot_grad(i) / radius - w(i));
    lambda0 = temp > lambda0 ? temp : lambda0; 
  }
  El::BigFloat lambda(lambda0);
  //cout << "1d search "<< endl;
  newton_search1d (thresh, max_it, w, lambda, e_dot_grad, radius, n_external_parameters);
  El::Zero(step);
  //El::BigFloat e_dot_var( - dot(Q(i, El::ALL),extvar));
  //cout << "axpy" <<endl;
  for (int i = 0; i < n_external_parameters ; i++){
     El::Axpy(( - El::Dot(Q(El::ALL,i),grad))/(w(i) + lambda), Q(El::ALL,i),step); 
  }  
}

void find_zero_step(const El::BigFloat &thresh, const int &max_it, const El::BigFloat &radius,
                                 El::Matrix<El::BigFloat> &hess, const El::Matrix<El::BigFloat> &grad, bool &find_zeros, 
                                 El::Matrix<El::BigFloat> &external_step, const El::BigFloat &quadratic_const){
  El::BigFloat quadratic_model = 0.5 * El::Dot(external_step, grad) + quadratic_const ;
  if (quadratic_model > 0) {
    find_zeros = false; 
    //cout << "hessian" << hess(0,0) << endl;
  }
  else{
    if (El::Nrm2(external_step) <= radius){
      find_zeros = true;
      //cout << "Looking for zeros using the original hessian" << endl;
    }
    else{
      trust_region_boundary_step (thresh, max_it, radius, 
                                  hess, grad, external_step);
      quadratic_model = 0.5 * El::Dot(external_step, grad) + quadratic_const;
      if (quadratic_model > 0) {
        find_zeros = false;
        cout << "hessian after correction: " << hess(0,0) << endl;
      }
      else{
        find_zeros = true;
        cout << "Looking for zeros using the updated hessian" << endl;
      }
    }
  }
}     
//
//int main(){
//  El::Environment env();
//  El::mpi::Comm comm = El::mpi::COMM_WORLD;
//  El::PushBlocksizeStack( 128 );
//  El::Matrix<El::BigFloat> testM;
//  El::Identity(testM,3,3);
//  //El::Jordan(testM,3,El::BigFloat(2));
//  //El::MakeHermitian(El::UPPER,testM);
//  testM(0,0) = 2;
//  testM(1,1) = 1;
//  testM(2,2) = 5;
//  testM(0,1) = 3; 
//  testM(0,2) = 0;
//  testM(1,2) = 1;
//  El::MakeHermitian(El::UPPER,testM);
//  //cout << "height: " << testM.Height() << endl;
//  El::Print(testM, "test matrix:");
//  El::Matrix<El::BigFloat> testG(3,1);
//  El::Zero(testG);
//  testG(0) = 0.3;
//  testG(1) = 0.2;
//  testG(2) = 0.4;
//  //El::Print(testG, "test grad:");
//  El::Matrix<El::BigFloat> testvar(3,1);
//  testvar(0) = 1.2;
//  testvar(1) = 0.1;
//  testvar(2) = 3;
//  //El::Print(testvar, "test var:");
//  El::Matrix<El::BigFloat> step(3,1); 
//  El::Zero(step);
//  //El::Print(step, "step:");
//  trust_region_boundary_step (0.0001,10,0.00001,testM, testG, testvar, step, 3);
//  //El::Print(step, "step:");
//  //      
//
//  return 0; 
//}

