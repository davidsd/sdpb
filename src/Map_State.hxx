#pragma once

#include <map>
#include <vector>
#include <string>
#include <stdexcept>

template <typename Key_State, typename Value_State> class Map_State
{
public:
  std::string name;
  bool inside = false, parsing_key = false, parsing_value = false,
       parsing_entry = false;
  ;
  Key_State key_state;
  Value_State value_state;
  std::map<decltype(key_state.value), decltype(value_state.value)> value;

  Map_State(const std::vector<std::string> &names, const size_t &offset)
      : name(names.at(offset)), key_state(names, offset + 1),
        value_state(names, offset + 2)
  {}
  Map_State(const std::initializer_list<std::string> &names)
      : Map_State(names, 0)
  {}

  Map_State() = delete;

  // JSON Functions
  void json_key(const std::string &key)
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_key(key);
              }
            else if(parsing_value)
              {
                value_state.json_key(key);
              }
            else
              {
                throw std::runtime_error(
                  "Invalid input file.  Found the key '" + key
                  + "' while parsing a key-value entry.");
              }
          }
        else
          {
            throw std::runtime_error("Invalid input file.  Found the key '"
                                     + key
                                     + "' while parsing a key-value map.");
          }
      }
    else
      {
        throw std::runtime_error(
          "Invalid input file.  Found the unexpected key '" + key + "'.");
      }
  }

  void json_string(const std::string &s)
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_string(s);
                parsing_key = key_state.inside;
              }
            else if(parsing_value)
              {
                value_state.json_string(s);
                parsing_value = value_state.inside;
              }
            else
              {
                throw std::runtime_error(
                  "Invalid input file.  Found the string '" + s
                  + "' while parsing a key-value entry.");
              }
          }
        else
          {
            throw std::runtime_error("Invalid input file.  Found the string '"
                                     + s + "' while parsing a key-value map.");
          }
      }
    else
      {
        throw std::runtime_error(
          "Invalid input file.  Found the unexpected string '" + s + "'.");
      }
  }

  void json_start_array()
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_start_array();
              }
            else if(parsing_value)
              {
                value_state.json_start_array();
              }
            else
              {
                throw std::runtime_error("Invalid input file.  Found an array "
                                         "inside a key-value entry.");
              }
          }
        else
          {
            parsing_entry = true;
            parsing_key=true;
          }
      }
    else
      {
        inside = true;
        value.clear();
      }
  }

  void json_end_array()
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_end_array();
                parsing_key=key_state.inside;
              }
            else if(parsing_value)
              {
                value_state.json_end_array();
                parsing_value=value_state.inside;
              }
            else
              {
                parsing_entry = false;
                value.emplace(key_state.value, value_state.value);
              }
          }
        else
          {
            inside = false;
          }
      }
    else
      {
        throw std::runtime_error(
          "Invalid input file.  Found an unexpected array.");
      }
  }

  void json_start_object()
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_start_object();
              }
            else if(parsing_value)
              {
                value_state.json_start_object();
              }
            else
              {
                throw std::runtime_error(
                  "Invalid input file.  Found an unexpected object start when "
                  "parsing a key-value entry.");
              }
          }
        else
          {
            throw std::runtime_error(
              "Invalid input file.  Found an unexpected object start when "
              "parsing a key-value map.");
          }
      }
    else
      {
        throw std::runtime_error(
          "Invalid input file.  Found an unexpected object start.");
      }
  }

  void json_end_object()
  {
    if(inside)
      {
        if(parsing_entry)
          {
            if(parsing_key)
              {
                key_state.json_end_object();
                parsing_key=key_state.inside;
              }
            else if(parsing_value)
              {
                value_state.json_end_object();
                parsing_value=value_state.inside;
              }
            else
              {
                throw std::runtime_error(
                  "Invalid input file.  Found an unexpected object end when "
                  "parsing a key-value entry.");
              }
          }
        else
          {
            throw std::runtime_error(
              "Invalid input file.  Found an unexpected object end when "
              "parsing a key-value map.");
          }
      }
    else
      {
        throw std::runtime_error(
          "Invalid input file.  Found an unexpected object end.");
      }
  }
};
