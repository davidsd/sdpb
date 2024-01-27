#pragma once

#include "Json_Polynomial_Vector_Parser.hxx"
#include "sdpb_util/json/Json_Vector_Parser.hxx"

using Vector_Of_Polynomial_Vectors = std::vector<Polynomial_Vector>;
using Matrix_Of_Polynomial_Vectors = std::vector<Vector_Of_Polynomial_Vectors>;

using Json_Vector_Of_Polynomial_Vectors_Parser
  = Json_Vector_Parser<Json_Polynomial_Vector_Parser>;

// Matrix = Vector of Vectors <Polynomial_Vector>
using Json_Matrix_Of_Polynomial_Vectors_Parser
  = Json_Vector_Parser<Json_Vector_Of_Polynomial_Vectors_Parser>;
