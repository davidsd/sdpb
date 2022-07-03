#include <El.hpp>


void print_matrix(El::Matrix<El::BigFloat> & matrix);

void BFGS_update_hessian(const int n_parameters, 
                           const El::Matrix<El::BigFloat> &grad_p_diff, 
                           const El::Matrix<El::BigFloat> &last_it_step, 
                           El::Matrix<El::BigFloat> &hess_bfgs
                           )
{
    El::Matrix<El::BigFloat> tempBs(last_it_step);
    El::Gemv(El::NORMAL, El::BigFloat(1), hess_bfgs, last_it_step,  El::BigFloat(0), tempBs);

    El::Matrix<El::BigFloat> tempyy(n_parameters,n_parameters);
    El::Zero(tempyy);
    El::Ger(El::BigFloat(1), grad_p_diff, grad_p_diff, tempyy);
    tempyy *= El::BigFloat(1)/El::Dot(grad_p_diff, last_it_step);

//	std::cout << "\ntempyy :\n";
//	print_matrix(tempyy);

    El::Matrix<El::BigFloat> tempBssB(n_parameters,n_parameters); 
    El::Zero(tempBssB);
    El::Ger(El::BigFloat(1),tempBs, tempBs, tempBssB);
    tempBssB *= El::BigFloat(1)/El::Dot(last_it_step, tempBs);

//	std::cout << "\ntempBssB :\n";
//	print_matrix(tempBssB);
    
//	std::cout << "hess_bfgs.Precision = " << hess_bfgs(0, 0).Precision();


    hess_bfgs += tempyy; 
    hess_bfgs -= tempBssB; 
}

//void BFGS_line_search(const El::Matrix<El::BigFloat> &search_direction, El::BigFloat scaling_factor)
//{
//     f(
//}
