#pragma once

#include "Damped_Rational_State.hxx"
#include "Polynomial_State.hxx"
#include "../Positive_Matrix_With_Prefactor.hxx"

#include <libxml2/libxml/parser.h>
#include <vector>
#include <string>
#include <stdexcept>

struct Positive_Matrix_With_Prefactor_State
{
  std::string name;
  bool inside = false, finished_damped_rational = false;
  Positive_Matrix_With_Prefactor value;

  Damped_Rational_State damped_rational_state;
  Vector_State<Vector_State<Vector_State<Polynomial_State>>> polynomials_state;

  Positive_Matrix_With_Prefactor_State(const std::vector<std::string> &names,
                                       const size_t &offset)
      : name(names.at(offset)), damped_rational_state(names, offset + 1),
        polynomials_state(names, offset + 2)
  {}
  Positive_Matrix_With_Prefactor_State(
    const std::initializer_list<std::string> &names)
      : Positive_Matrix_With_Prefactor_State(names, 0)
  {}

  bool on_start_element(const std::string &element_name)
  {
    if(inside)
      {
        if(element_name != "Symbol"
           && (finished_damped_rational
               || !damped_rational_state.on_start_element(element_name))
           && !polynomials_state.on_start_element(element_name))
          {
            throw std::runtime_error(
              "Invalid input file.  Expected 'Function' inside "
              "PositiveMatrixWithPrefactor, but found '"
              + element_name + "'");
          }
      }
    else
      {
        inside = (element_name == name);
        if(inside)
          {
            finished_damped_rational = false;
          }
      }
    return inside;
  }

  bool on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(damped_rational_state.on_end_element(element_name))
          {
            finished_damped_rational = !damped_rational_state.inside;
            if(finished_damped_rational)
              {
                using namespace std;
                swap(damped_rational_state.value, value.damped_rational);
              }
          }
        else if(polynomials_state.on_end_element(element_name))
          {
            if(!polynomials_state.inside)
              {
                value.polynomials.reserve(polynomials_state.value.size());
                for(auto &polynomial_vector_vector : polynomials_state.value)
                  {
                    value.polynomials.emplace_back();
                    auto &vector_vector_poly(value.polynomials.back());
                    vector_vector_poly.reserve(polynomial_vector_vector.size());
                    for(auto &polynomial_vector: polynomial_vector_vector)
                      {
                        vector_vector_poly.emplace_back();
                        auto &vector_poly(vector_vector_poly.back());
                        vector_poly.reserve(polynomial_vector.size());
                        for(auto &polynomial: polynomial_vector)
                          {
                            vector_poly.emplace_back();
                            auto &poly(vector_poly.back());
                            poly.coefficients.resize(polynomial.size(),0);
                            for(auto &term : polynomial)
                              {
                                if(poly.coefficients.size() < term.first + 1)
                                  {
                                    poly.coefficients.resize(term.first + 1, 0);
                                  }
                                poly.coefficients[term.first] = term.second;
                              }
                          }
                      }
                  }
              }
          }
        else
          {
            inside = (element_name != name);
          }
      }
    return result;
  }

  bool on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        damped_rational_state.on_characters(characters, length)
          || polynomials_state.on_characters(characters, length);
      }
    return inside;
  }
};
