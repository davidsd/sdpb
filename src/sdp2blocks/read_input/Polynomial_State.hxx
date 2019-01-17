#pragma once

#include "Polynomial_Term_State.hxx"

struct Polynomial_State
{
  std::string name;
  bool inside = false;
  std::vector<std::pair<size_t, El::BigFloat>> value;

  Number_State number_state;
  Vector_State<Polynomial_Term_State> vector_polynomial_term_state;

  Polynomial_State(const std::vector<std::string> &names, const size_t &offset)
      : name(names.at(offset)), number_state("Number"s),
        vector_polynomial_term_state({"Function"s,"Function"s,})
  {}
  Polynomial_State(const std::initializer_list<std::string> &names)
      : Polynomial_State(names, 0)
  {}

  bool on_start_element(const std::string &element_name)
  {
    // if(element_name=="Number")
    //   std::cout << "Polynomial Starting\n";
    
    if((vector_polynomial_term_state.inside
        || !number_state.on_start_element(element_name))
       && !vector_polynomial_term_state.on_start_element(element_name))
      {
        throw std::runtime_error(
          "Invalid input file.  Expected '" + vector_polynomial_term_state.name
          + "' or '" + number_state.name + "' inside Polynomial, but found '"
          + element_name + "'");
      }
    inside = (number_state.inside || vector_polynomial_term_state.inside);

    // if(element_name=="Number")
    //   std::cout << "Polynomial Started\n";
    return inside;
  }

  bool on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(number_state.on_end_element(element_name))
          {
            if(!number_state.inside)
              {
                value.emplace_back(0, number_state.value);
              }
          }
        else if(vector_polynomial_term_state.on_end_element(element_name))
          {
            if(!vector_polynomial_term_state.inside)
              {
                using namespace std;
                swap(value, vector_polynomial_term_state.value);
              }
          }
        inside = (number_state.inside || vector_polynomial_term_state.inside);
      }
    return result;
  }

  bool on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        number_state.on_characters(characters, length)
          || vector_polynomial_term_state.on_characters(characters, length);
      }
    return inside;
  }
};
