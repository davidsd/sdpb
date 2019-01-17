#pragma once

#include "Damped_Rational_State.hxx"
#include "Polynomial_Term_State.hxx"
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
  using Polynomial_State = Vector_State<Polynomial_Term_State>;
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
      }
    // std::cout << "Positive start: " << inside << " "
    //           << damped_rational_state.inside << " "
    //           << polynomials_state.inside << "\n";
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
                // std::cout
                //   << "sizes: "
                //   << polynomials_state.value.front().front().front().size()
                //   << " " << polynomials_state.value.front().front().size()
                //   << " " << polynomials_state.value.front().size() << " "
                //   << polynomials_state.value.size() << " "
                //   << "\n";
                for(auto &polynomial : polynomials_state.value.front().front())
                  {
                    std::vector<El::BigFloat> coefficients(polynomial.size(),
                                                           0);
                    for(auto &term : polynomial)
                      {
                        if(coefficients.size() < term.first + 1)
                          {
                            coefficients.resize(term.first + 1, 0);
                          }
                        coefficients[term.first] = term.second;
                      }
                    value.polynomials.emplace_back();
                    swap(value.polynomials.back().coefficients, coefficients);
                  }
              }
          }
        else
          {
            inside = (element_name != name);
          }
      }
    // std::cout << "Positive end: " << inside << " "
    //           << damped_rational_state.inside << " "
    //           << polynomials_state.inside << "\n";
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
