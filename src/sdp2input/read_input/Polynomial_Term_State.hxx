#pragma once

#include "Times_State.hxx"

struct Polynomial_Term_State
{
  std::string name;
  bool inside = false;
  std::pair<size_t, El::BigFloat> value;

  Number_State<El::BigFloat> number_state;
  Times_State times_state;

  Polynomial_Term_State(const std::vector<std::string> &names,
                        const size_t &offset)
      : name(names.at(offset)), number_state("Number"), times_state("Function")
  {}

  bool xml_on_start_element(const std::string &element_name)
  {
    // Have to check times_state before number_state, because inside
    // times_state there is a valid Number.
    if(element_name != "Symbol"
       && (times_state.inside
           || !number_state.xml_on_start_element(element_name))
       && !times_state.xml_on_start_element(element_name))
      {
        throw std::runtime_error(
          "Invalid input file.  Expected 'Number' or 'Function' inside "
          "Polynomial_Term, but found '"
          + element_name + "'");
      }
    inside = (number_state.inside || times_state.inside);
    return inside;
  }

  bool xml_on_end_element(const std::string &element_name)
  {
    bool result(inside);
    if(inside)
      {
        if(number_state.xml_on_end_element(element_name))
          {
            if(!number_state.inside)
              {
                value.first = 0;
                using namespace std;
                swap(number_state.value, value.second);
              }
          }
        else if(times_state.xml_on_end_element(element_name))
          {
            if(!times_state.inside)
              {
                using namespace std;
                swap(value, times_state.value);
              }
          }
        inside = (number_state.inside || times_state.inside);
      }
    return result;
  }

  bool xml_on_characters(const xmlChar *characters, int length)
  {
    if(inside)
      {
        number_state.xml_on_characters(characters, length)
          || times_state.xml_on_characters(characters, length);
      }
    return inside;
  }
};
