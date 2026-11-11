// GPU backend for SDPB: a device float on the host, the value that crosses the
// boundary between the device code and SDPB's El::BigFloat driver.
#pragma once

#include <cstdint>
#include <vector>

namespace sdpb_gpu
{
  // value = m * 2^e, m two's complement, 32-bit limbs, little endian; zero:
  // no limbs or all limbs zero
  struct Host_Float
  {
    std::vector<uint32_t> m;
    int64_t e = -(int64_t(1) << 40);
  };
}
