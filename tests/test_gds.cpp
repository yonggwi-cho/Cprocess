#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "cprocess/gds_reader.hpp"
#include "test_util.hpp"

using namespace cp;

namespace {

// ---------------------------------------------------------------------
// Minimal test-local GDSII byte-level writer (not the production code).
// ---------------------------------------------------------------------

void put_u16(std::vector<unsigned char>& b, uint16_t v) {
  b.push_back(static_cast<unsigned char>(v >> 8));
  b.push_back(static_cast<unsigned char>(v & 0xFF));
}
void put_i32(std::vector<unsigned char>& b, int32_t v) {
  uint32_t u = static_cast<uint32_t>(v);
  b.push_back(static_cast<unsigned char>((u >> 24) & 0xFF));
  b.push_back(static_cast<unsigned char>((u >> 16) & 0xFF));
  b.push_back(static_cast<unsigned char>((u >> 8) & 0xFF));
  b.push_back(static_cast<unsigned char>(u & 0xFF));
}

// rectype/datatype combined per the spec's hex record codes (e.g. 0x0102 ->
// rectype 0x01, datatype 0x02).
void record(std::vector<unsigned char>& b, uint16_t code,
            const std::vector<unsigned char>& payload) {
  uint16_t len = static_cast<uint16_t>(4 + payload.size());
  put_u16(b, len);
  b.push_back(static_cast<unsigned char>(code >> 8));
  b.push_back(static_cast<unsigned char>(code & 0xFF));
  b.insert(b.end(), payload.begin(), payload.end());
}

void record_i16(std::vector<unsigned char>& b, uint16_t code,
                const std::vector<int16_t>& vals) {
  std::vector<unsigned char> p;
  for (int16_t v : vals) put_u16(p, static_cast<uint16_t>(v));
  record(b, code, p);
}

void record_i32(std::vector<unsigned char>& b, uint16_t code,
                const std::vector<int32_t>& vals) {
  std::vector<unsigned char> p;
  for (int32_t v : vals) put_i32(p, v);
  record(b, code, p);
}

void record_str(std::vector<unsigned char>& b, uint16_t code,
                const std::string& s) {
  std::string padded = s;
  if (padded.size() % 2 != 0) padded.push_back('\0');
  std::vector<unsigned char> p(padded.begin(), padded.end());
  record(b, code, p);
}

void record_empty(std::vector<unsigned char>& b, uint16_t code) {
  record(b, code, {});
}

// UNITS payload: two 8-byte excess-64 reals, given verbatim per the spec
// (1e-3 user-units-per-DB-unit, 1e-9 meters-per-DB-unit).
void record_units(std::vector<unsigned char>& b) {
  std::vector<unsigned char> p = {
      0x3E, 0x41, 0x89, 0x37, 0x4B, 0xC6, 0xA7, 0xF0,  // 1e-3
      0x39, 0x44, 0xB8, 0x2F, 0xA0, 0x9B, 0x5A, 0x54,  // 1e-9
  };
  record(b, 0x0305, p);
}

std::string write_tmp(const char* name, const std::vector<unsigned char>& bytes) {
  std::string path = std::string("./") + name;
  std::ofstream f(path, std::ios::binary);
  f.write(reinterpret_cast<const char*>(bytes.data()),
          static_cast<std::streamsize>(bytes.size()));
  return path;
}

// Builds a GDS stream with:
//  - BOUNDARY on layer 2: 1um square (0,0)-(1000,0)-(1000,1000)-(0,1000),
//    closed with a duplicate first point.
//  - a PATH record (unrecognized -> must be skipped without crashing)
//  - BOUNDARY on layer 5: triangle (0,0)-(500,0)-(250,500), closed.
std::vector<unsigned char> build_test_gds() {
  std::vector<unsigned char> b;
  record_empty(b, 0x0102);              // HEADER (payload omitted for brevity)
  record_i16(b, 0x0202, std::vector<int16_t>(12, 0));  // BGNLIB
  record_str(b, 0x0206, "LIB");         // LIBNAME
  record_units(b);                      // UNITS
  record_i16(b, 0x0502, std::vector<int16_t>(12, 0));  // BGNSTR
  record_str(b, 0x0606, "TOP");         // STRNAME

  // BOUNDARY layer 2: square
  record_empty(b, 0x0800);              // BOUNDARY
  record_i16(b, 0x0D02, {2});           // LAYER
  record_i16(b, 0x0E02, {0});           // DATATYPE
  record_i32(b, 0x1003,
             {0, 0, 1000, 0, 1000, 1000, 0, 1000, 0, 0});  // XY (closed)
  record_empty(b, 0x1100);              // ENDEL

  // A PATH record (unrecognized type here) mixed in -- must be skipped.
  {
    std::vector<unsigned char> p;
    put_i32(p, 0);
    put_i32(p, 0);
    put_i32(p, 100);
    put_i32(p, 100);
    record(b, 0x0906, p);  // pretend PATH-ish record, not in our switch table
  }

  // BOUNDARY layer 5: triangle
  record_empty(b, 0x0800);              // BOUNDARY
  record_i16(b, 0x0D02, {5});           // LAYER
  record_i16(b, 0x0E02, {0});           // DATATYPE
  record_i32(b, 0x1003, {0, 0, 500, 0, 250, 500, 0, 0});  // XY (closed)
  record_empty(b, 0x1100);              // ENDEL

  record_empty(b, 0x0700);              // ENDSTR
  record_empty(b, 0x0400);              // ENDLIB
  return b;
}

}  // namespace

int main() {
  // 1. Roundtrip + 2. layer filter + 4. closed-polygon dedup.
  {
    auto bytes = build_test_gds();
    std::string path = write_tmp("t_gds_ok.gds", bytes);

    auto sq = read_gds(path, 2);
    CHECK(sq.size() == 1);
    CHECK(sq[0].size() == 4);  // dedup: no trailing duplicate point
    CHECK_NEAR(sq[0][0].first, 0.0, 1e-9);
    CHECK_NEAR(sq[0][0].second, 0.0, 1e-9);
    CHECK_NEAR(sq[0][1].first, 1e-4, 1e-9);   // 1000 DBU * 1e-9 m/DBU * 100 = 1e-4 cm
    CHECK_NEAR(sq[0][1].second, 0.0, 1e-9);
    CHECK_NEAR(sq[0][2].first, 1e-4, 1e-9);
    CHECK_NEAR(sq[0][2].second, 1e-4, 1e-9);
    CHECK_NEAR(sq[0][3].first, 0.0, 1e-9);
    CHECK_NEAR(sq[0][3].second, 1e-4, 1e-9);

    auto tri = read_gds(path, 5);
    CHECK(tri.size() == 1);
    CHECK(tri[0].size() == 3);  // dedup applied here too

    auto none = read_gds(path, 1);
    CHECK(none.empty());

    auto all = read_gds(path, -1);
    CHECK(all.size() == 2);

    std::remove(path.c_str());
  }

  // 3. Error cases.
  {
    bool threw = false;
    try {
      read_gds("./no_such_file_xyz123.gds", -1);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
  }
  {
    std::vector<unsigned char> junk = {'n', 'o', 't', ' ', 'a', ' ', 'g', 'd',
                                        's', ' ', 'f', 'i', 'l', 'e', '\n'};
    std::string path = write_tmp("t_gds_bad.gds", junk);
    bool threw = false;
    try {
      read_gds(path, -1);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    CHECK(threw);
    std::remove(path.c_str());
  }

  std::printf("gds tests passed\n");
  return 0;
}
