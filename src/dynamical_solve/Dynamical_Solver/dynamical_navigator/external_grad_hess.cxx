#include <El.hpp>

//Given objectives on a grid in the external parameter (p) space,
//where ePlus correspond to (+e_i), eMinus corresponds to (-e_i) and eSum corresponds to (e_i+e_j) directions respectively, 
//Compute H_pp and Del_p(L_mu) to solve Eq(13).
//Return: void.
//Update: hess = H_pp ,
//        grad = Del_p(L_mu).
void external_grad_hessian(const El::Matrix<El::BigFloat> &ePlus, 
                           const El::Matrix<El::BigFloat> &eMinus, 
                           const El::Matrix<El::BigFloat> &eSum,
                           const El::BigFloat &alpha,
                           El::Matrix<El::BigFloat> &grad,
                           El::Matrix<El::BigFloat> &hess)
{
  grad =  ePlus;
  grad -= eMinus;
  grad *= 1.0/(2.0*alpha);
  hess = eSum;
  for (int i=0; i<ePlus.Height(); i++)
    { 
      hess(i,i) += (ePlus(i) + eMinus(i) );
      for (int j=0; j<i ; j++)
        {
          hess(i,j) += (- ePlus(i) - ePlus(j));
          hess(j,i) += (- ePlus(i) - ePlus(j));
        }
    }
  hess *= 1.0/(alpha * alpha);
}