#include "integration_tests/common.hxx"
#include "sdpa_solve/parse_sdpa/parse_sdpa.hxx"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
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

        // The visitor runs on several threads; collect per slot, then apply
        // the entries in file order.
        struct Seen
        {
          uint64_t order;
          Key key;
          El::BigFloat value;
        };
        std::vector<std::vector<Seen>> seen(sdpa_visit_slots());
        visit_sdpa_lower(file.path(), [&](size_t slot, uint64_t order,
                                         size_t block, size_t matrix, size_t i,
                                         size_t j, const El::BigFloat &value) {
          seen.at(slot).push_back({order, Key{block, matrix, i, j}, value});
        });
        std::vector<Seen> all;
        for(auto &slot : seen)
          for(auto &entry : slot)
            all.push_back(std::move(entry));
        std::sort(all.begin(), all.end(), [](const Seen &a, const Seen &b) {
          return a.order < b.order;
        });
        for(size_t k = 1; k < all.size(); ++k)
          REQUIRE(all[k - 1].order < all[k].order);
        for(const auto &entry : all)
          {
            const auto [block, matrix, i, j] = entry.key;
            REQUIRE(i >= j);
            if(entry.value != zero)
              actual[entry.key] = entry.value;
            else
              actual.erase(entry.key);
          }
        REQUIRE(!expected.empty());
        REQUIRE(actual == expected);
      }
}

// A sparse file of several reader slices (16 MB each), read in more than one
// batch of slices: many repeated positions, explicit zeros that cancel earlier
// values, and lines cut by slice ends. Plain and gzip input.
TEST_CASE("sdpa_streamed_input_slices")
{
  const auto dir = Test_Config::test_output_dir() / "sdpa_streamed_input_slices";
  std::filesystem::create_directories(dir);
  const std::vector<int> sizes{12, -4, 7};
  const size_t m = 9;
  std::string text = "\"sdpb.objectiveConst 1.5\n" + std::to_string(m) + "\n3\n12 -4 7\n";
  for(size_t p = 1; p <= m; ++p)
    text += std::to_string(p) + " ";
  text += "\n";
  uint64_t state = 12345;
  auto next = [&](uint64_t bound) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (state >> 33) % bound;
  };
  while(text.size() < (size_t(40) << 20))
    {
      const size_t p = next(m + 1), b = next(3);
      const size_t n = std::abs(sizes[b]);
      size_t i = 1 + next(n), j = sizes[b] < 0 ? i : 1 + next(n);
      if(i > j)
        std::swap(i, j);
      std::string value = next(5) == 0 ? "0"
                          : (next(2) ? "-" : "") + std::to_string(next(1000000))
                              + "." + std::to_string(next(1000000)) + "e"
                              + std::to_string(int(next(40)) - 20);
      text += std::to_string(p) + " " + std::to_string(b + 1) + " "
              + std::to_string(i) + " " + std::to_string(j) + " " + value + "\n";
    }
  const auto plain = dir / "slices.dat-s", gzip = dir / "slices.dat-s.gz";
  std::ofstream(plain, std::ios::binary) << text;
  {
    std::ofstream out(gzip, std::ios::binary);
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_compressor());
    os.push(out);
    os << text;
  }

  const auto cpu = read_sdpa(plain, [](size_t) { return true; });
  const El::BigFloat zero(0);
  using Key = std::tuple<size_t, size_t, size_t, size_t>;
  std::map<Key, El::BigFloat> expected;
  for(const auto &[block, matrices] : cpu.parsed_blocks)
    for(size_t p = 0; p < matrices.size(); ++p)
      for(int i = 0; i < matrices[p].Height(); ++i)
        for(int j = 0; j <= i; ++j)
          if(matrices[p](i, j) != zero)
            expected[{block, p, i, j}] = matrices[p](i, j);
  REQUIRE(!expected.empty());

  const char *threads = std::getenv("OMP_NUM_THREADS");
  const std::string saved = threads ? threads : "";
  for(const char *nt : {"1", "2", "5"})
    for(const auto &file : {plain, gzip})
      DYNAMIC_SECTION(file.filename().string() << ", " << nt << " threads")
      {
        setenv("OMP_NUM_THREADS", nt, 1);
        std::vector<std::vector<std::pair<uint64_t, std::pair<Key, El::BigFloat>>>>
          seen(sdpa_visit_slots());
        visit_sdpa_lower(file, [&](size_t slot, uint64_t order, size_t block,
                                  size_t matrix, size_t i, size_t j,
                                  const El::BigFloat &value) {
          seen.at(slot).push_back({order, {Key{block, matrix, i, j}, value}});
        });
        if(threads)
          setenv("OMP_NUM_THREADS", saved.c_str(), 1);
        else
          unsetenv("OMP_NUM_THREADS");
        std::vector<std::pair<uint64_t, std::pair<Key, El::BigFloat>>> all;
        size_t used_slots = 0;
        for(auto &slot : seen)
          {
            used_slots += !slot.empty();
            for(auto &entry : slot)
              all.push_back(std::move(entry));
          }
        // About 40 MB make three slices
        REQUIRE(used_slots == std::min<size_t>(seen.size(), 3));
        std::sort(all.begin(), all.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });
        size_t lines = std::count(text.begin(), text.end(), '\n') - 5;
        REQUIRE(all.size() == lines);
        std::map<Key, El::BigFloat> actual;
        for(size_t k = 0; k < all.size(); ++k)
          {
            REQUIRE((k == 0 || all[k - 1].first < all[k].first));
            const auto &[key, value] = all[k].second;
            if(value != zero)
              actual[key] = value;
            else
              actual.erase(key);
          }
        REQUIRE(actual == expected);
      }
  std::filesystem::remove_all(dir);
}
