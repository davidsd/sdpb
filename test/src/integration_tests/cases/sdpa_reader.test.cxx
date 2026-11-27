#include "integration_tests/common.hxx"
#include "sdpa_solve/parse_sdpa/parse_sdpa.hxx"

#include <map>
#include <tuple>

using namespace Sdpb::Sdpa;
using namespace Test_Util;

// Compare the sparse input representation used by the GPU with the full
// matrices used by the CPU, including implicit zeros and split diagonal blocks.
TEST_CASE("sdpa_streamed_input")
{
  const auto fixtures = Test_Config::test_data_dir / "end-to-end_tests" / "sdpa";
  for(const auto &problem : std::filesystem::directory_iterator(fixtures))
    for(const auto &file : std::filesystem::directory_iterator(problem.path() / "input"))
      DYNAMIC_SECTION(file.path().string())
      {
        const auto cpu = read_sdpa(file.path(), [](size_t) { return true; });
        const auto objective = read_sdpa_objective(file.path());
        REQUIRE(objective.c_objective == cpu.c_objective);
        REQUIRE(objective.objective_const == cpu.objective_const);

        const El::BigFloat zero(0);
        using Key = std::tuple<size_t, size_t, size_t, size_t>;
        std::map<Key, El::BigFloat> expected, actual;
        for(const auto &[block, matrices] : cpu.parsed_blocks)
          for(size_t p = 0; p < matrices.size(); ++p)
            for(int i = 0; i < matrices[p].Height(); ++i)
              for(int j = 0; j <= i; ++j)
                if(matrices[p](i, j) != zero)
                  expected[{block, p, i, j}] = matrices[p](i, j);

        visit_sdpa_lower(file.path(), [&](size_t block, size_t matrix, size_t i,
                                         size_t j, const El::BigFloat &value) {
          REQUIRE(i >= j);
          const Key key{block, matrix, i, j};
          if(value != zero)
            actual[key] = value;
          else
            actual.erase(key);
        });
        REQUIRE(!expected.empty());
        REQUIRE(actual == expected);
      }
}
