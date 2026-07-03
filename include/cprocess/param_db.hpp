#pragma once
#include <map>
#include <string>

namespace cp {

// Runtime parameter-override registry (singleton). Lets materials.cpp (and a
// few TED constants in diffusion.cpp/process.cpp) be tuned at runtime without
// rebuilding, for calibration and sensitivity analysis.
//
// Thread safety: set()/erase()/clear() must NOT be called while a solve is in
// progress (no lock is taken). Control flow is single-threaded; the OpenMP
// parallel regions inside solve() only call get() (read-only), which is safe
// concurrently since std::map lookups don't mutate the container.
//
// Key naming: "<Sym>.<field>" for per-dopant Arrhenius parameters (e.g.
// "B.d0", "B.seg_m0"), "I.<field>" for the self-interstitial point-defect
// model, "ted.<field>" for TED-loop constants. See the key table in
// materials.hpp for the full supported set. Unknown keys are accepted
// without validation -- they are simply inert unless some call site reads
// them.
class ParamDB {
 public:
  static ParamDB& instance();

  // Returns the override for `key` if set, else `fallback`.
  double get(const std::string& key, double fallback) const;
  void set(const std::string& key, double value);
  bool erase(const std::string& key);  // true if a value was removed
  void clear();                        // remove all overrides (tests)
  std::map<std::string, double> all() const;

 private:
  ParamDB() = default;
  std::map<std::string, double> overrides_;
};

}  // namespace cp
