#include "cprocess/param_db.hpp"

namespace cp {

ParamDB& ParamDB::instance() {
  static ParamDB db;
  return db;
}

double ParamDB::get(const std::string& key, double fallback) const {
  auto it = overrides_.find(key);
  return (it != overrides_.end()) ? it->second : fallback;
}

void ParamDB::set(const std::string& key, double value) {
  overrides_[key] = value;
}

bool ParamDB::erase(const std::string& key) {
  return overrides_.erase(key) > 0;
}

void ParamDB::clear() { overrides_.clear(); }

std::map<std::string, double> ParamDB::all() const { return overrides_; }

}  // namespace cp
