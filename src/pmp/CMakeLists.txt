add_library(pmp SHARED)
target_sources(pmp
                PRIVATE
                Polynomial_Vector_Matrix.cxx
                convert/sample_points.cxx
                convert/sample_scalings.cxx
                convert/bilinear_basis/bilinear_basis.cxx
                convert/bilinear_basis/precompute/precompute.cxx
                convert/bilinear_basis/precompute/integral.cxx
                convert/bilinear_basis/bilinear_form/bilinear_form.cxx
                convert/bilinear_basis/bilinear_form/rest.cxx
                convert/bilinear_basis/bilinear_form/dExp.cxx
                convert/bilinear_basis/bilinear_form/derivative.cxx
                convert/bilinear_basis/bilinear_form/operator_plus_set_Derivative_Term.cxx)
target_link_libraries(pmp PUBLIC ${USE_PACKAGES})