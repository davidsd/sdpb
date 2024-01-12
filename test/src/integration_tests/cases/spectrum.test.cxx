#include "integration_tests/common.hxx"

#include <filesystem>

namespace fs = std::filesystem;

TEST_CASE("spectrum")
{
  INFO("Simple spectrum test");

  for(std::string name : {"1d", "matrix"})
    DYNAMIC_SECTION(name)
    {
      Test_Util::Test_Case_Runner runner("spectrum/" + name);

      fs::path data_dir = runner.data_dir;
      fs::path output_dir = runner.output_dir;

      Test_Util::Test_Case_Runner::Named_Args_Map args{
        {"--input", (data_dir / "pvm.xml").string()},
        {"--solution", (data_dir / "solution").string()},
        {"--output", (output_dir / "spectrum.json").string()},
        {"--precision", "1024"},
        {"--threshold", "1e-10"},
        // --format is obsolete, we keep to here to check backward compatibility,
        // i.e. spectrum shouldn't fail when we pass this option.
        {"--format", "PVM"}};

      runner.mpi_run({"build/spectrum"}, args, 2);

      auto out = output_dir / "spectrum.json";
      auto out_orig = data_dir / "spectrum.orig.json";
      Test_Util::REQUIRE_Equal::diff_spectrum(out, out_orig, 1024, 1024 / 2);
    }
}
