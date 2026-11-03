add_library(pmp STATIC
  Polynomial_Vector_Matrix.cxx
  convert/sample_points.cxx
  convert/sample_scalings.cxx
  convert/bilinear_basis/bilinear_basis.cxx)
target_link_libraries(pmp PUBLIC sdpb_util)
