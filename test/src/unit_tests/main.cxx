#include <catch2/catch_amalgamated.hpp>
#include <El.hpp>

#ifndef CATCH_AMALGAMATED_CUSTOM_MAIN
#error "To override main, pass '-D CATCH_AMALGAMATED_CUSTOM_MAIN' to compiler"
#endif

int main(int argc, char *argv[])
{
  // your setup ...
  El::Environment env(argc, argv);
  El::gmp::SetPrecision(128);

  int result = Catch::Session().run(argc, argv);

  // your clean-up...

  return result;
}
