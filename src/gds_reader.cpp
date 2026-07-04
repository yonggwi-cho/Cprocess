#include "cprocess/gds_reader.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace cp {
namespace {

constexpr uint16_t kHEADER = 0x0102;
constexpr uint16_t kBGNLIB = 0x0202;
constexpr uint16_t kLIBNAME = 0x0206;
constexpr uint16_t kUNITS = 0x0305;
constexpr uint16_t kBGNSTR = 0x0502;
constexpr uint16_t kSTRNAME = 0x0606;
constexpr uint16_t kENDSTR = 0x0700;
constexpr uint16_t kBOUNDARY = 0x0800;
constexpr uint16_t kLAYER = 0x0D02;
constexpr uint16_t kDATATYPE = 0x0E02;
constexpr uint16_t kXY = 0x1003;
constexpr uint16_t kENDEL = 0x1100;
constexpr uint16_t kENDLIB = 0x0400;

uint16_t u16be(const unsigned char* p) {
  return static_cast<uint16_t>((p[0] << 8) | p[1]);
}

int32_t i32be(const unsigned char* p) {
  uint32_t v = (static_cast<uint32_t>(p[0]) << 24) |
               (static_cast<uint32_t>(p[1]) << 16) |
               (static_cast<uint32_t>(p[2]) << 8) |
               static_cast<uint32_t>(p[3]);
  return static_cast<int32_t>(v);
}

int16_t i16be(const unsigned char* p) {
  return static_cast<int16_t>(u16be(p));
}

// GDSII 8-byte "excess-64" real: sign(1b) | exponent(7b, base-16, bias 64) |
// mantissa(56b). value = (-1)^s * mantissa/2^56 * 16^(exp-64).
double excess64_to_double(const unsigned char* p) {
  uint64_t raw = 0;
  for (int i = 0; i < 8; ++i) raw = (raw << 8) | p[i];
  bool sign = (raw >> 63) & 0x1;
  int exponent = static_cast<int>((raw >> 56) & 0x7F);
  uint64_t mantissa = raw & 0x00FFFFFFFFFFFFFFULL;
  if (mantissa == 0 && exponent == 0) return 0.0;
  double m = static_cast<double>(mantissa) / 72057594037927936.0;  // 2^56
  double value = m * std::pow(16.0, exponent - 64);
  return sign ? -value : value;
}

}  // namespace

std::vector<std::vector<std::pair<double, double>>> read_gds(
    const std::string& path, int layer) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("read_gds: cannot open file: " + path);

  std::vector<unsigned char> data((std::istreambuf_iterator<char>(in)),
                                   std::istreambuf_iterator<char>());
  if (data.size() < 4)
    throw std::runtime_error("read_gds: file too small / not a GDSII stream");

  double meters_per_dbu = 1e-9;  // fallback default
  bool have_header = false;
  bool in_boundary = false;
  int cur_layer = -1;
  std::vector<std::pair<double, double>> cur_xy;

  std::vector<std::vector<std::pair<double, double>>> result;

  size_t pos = 0;
  size_t n = data.size();
  while (pos < n) {
    if (pos + 4 > n)
      throw std::runtime_error("read_gds: truncated record header");
    uint16_t reclen = u16be(&data[pos]);
    if (reclen < 4)
      throw std::runtime_error("read_gds: invalid record length");
    if (pos + reclen > n)
      throw std::runtime_error("read_gds: truncated record payload");
    unsigned char rectype = data[pos + 2];
    unsigned char datatype = data[pos + 3];
    uint16_t code = static_cast<uint16_t>((rectype << 8) | datatype);
    const unsigned char* payload = &data[pos + 4];
    size_t paylen = reclen - 4;

    switch (code) {
      case kHEADER:
        have_header = true;
        break;
      case kBGNLIB:
      case kLIBNAME:
      case kBGNSTR:
      case kSTRNAME:
      case kENDSTR:
        break;
      case kUNITS: {
        if (paylen >= 16) {
          // 2nd value = meters per database unit.
          meters_per_dbu = excess64_to_double(payload + 8);
        }
        break;
      }
      case kBOUNDARY:
        in_boundary = true;
        cur_layer = -1;
        cur_xy.clear();
        break;
      case kLAYER:
        if (in_boundary && paylen >= 2) cur_layer = i16be(payload);
        break;
      case kDATATYPE:
        break;
      case kXY: {
        if (in_boundary) {
          size_t npts = paylen / 8;
          cur_xy.clear();
          cur_xy.reserve(npts);
          for (size_t i = 0; i < npts; ++i) {
            int32_t xi = i32be(payload + i * 8);
            int32_t yi = i32be(payload + i * 8 + 4);
            double x_cm = xi * meters_per_dbu * 100.0;
            double y_cm = yi * meters_per_dbu * 100.0;
            cur_xy.emplace_back(x_cm, y_cm);
          }
        }
        break;
      }
      case kENDEL: {
        if (in_boundary) {
          // Strip trailing duplicate closing point.
          if (cur_xy.size() >= 2 && cur_xy.front() == cur_xy.back())
            cur_xy.pop_back();
          if (layer < 0 || cur_layer == layer) {
            if (!cur_xy.empty()) result.push_back(cur_xy);
          }
        }
        in_boundary = false;
        cur_xy.clear();
        cur_layer = -1;
        break;
      }
      case kENDLIB:
        break;
      default:
        // SREF/AREF/PATH/TEXT/BOX/... and anything else: skip payload.
        break;
    }

    pos += reclen;
  }

  if (!have_header)
    throw std::runtime_error("read_gds: missing HEADER record (not a valid GDSII stream)");

  return result;
}

}  // namespace cp
