#include "cprocess/vtk_writer.hpp"

#include <cstdio>
#include <fstream>
#include <stdexcept>

namespace cp {

void write_vtu(
    const std::string& path, const Mesh& mesh,
    const std::vector<std::pair<std::string, const std::vector<double>*>>& scalars,
    const std::vector<std::pair<std::string, const std::vector<int>*>>& int_scalars,
    const std::vector<std::pair<std::string, const std::vector<double>*>>&
        point_scalars) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("vtu: cannot open for writing: " + path);

  const std::size_t np = mesh.nodes.size(), ncl = mesh.cells.size();
  const double to_um = 1e4;  // cm -> um for display
  char buf[128];

  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
         "byte_order=\"LittleEndian\">\n"
      << "<UnstructuredGrid>\n"
      << "<Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << ncl
      << "\">\n";

  out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";
  for (const auto& p : mesh.nodes) {
    std::snprintf(buf, sizeof(buf), "%.9g %.9g %.9g\n", p.x * to_um,
                  p.y * to_um, p.z * to_um);
    out << buf;
  }
  out << "</DataArray>\n</Points>\n";

  out << "<Cells>\n<DataArray type=\"Int64\" Name=\"connectivity\" "
         "format=\"ascii\">\n";
  for (const auto& c : mesh.cells)
    out << c[0] << ' ' << c[1] << ' ' << c[2] << ' ' << c[3] << '\n';
  out << "</DataArray>\n<DataArray type=\"Int64\" Name=\"offsets\" "
         "format=\"ascii\">\n";
  for (std::size_t i = 1; i <= ncl; ++i) out << 4 * i << '\n';
  out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" "
         "format=\"ascii\">\n";
  for (std::size_t i = 0; i < ncl; ++i) out << "10\n";  // VTK_TETRA
  out << "</DataArray>\n</Cells>\n";

  if (!point_scalars.empty()) {
    out << "<PointData>\n";
    for (const auto& [name, data] : point_scalars) {
      out << "<DataArray type=\"Float64\" Name=\"" << name
          << "\" format=\"ascii\">\n";
      for (std::size_t i = 0; i < np; ++i) {
        std::snprintf(buf, sizeof(buf), "%.6e\n",
                      i < data->size() ? (*data)[i] : 0.0);
        out << buf;
      }
      out << "</DataArray>\n";
    }
    out << "</PointData>\n";
  }

  out << "<CellData>\n";
  for (const auto& [name, data] : scalars) {
    out << "<DataArray type=\"Float64\" Name=\"" << name
        << "\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < ncl; ++i) {
      std::snprintf(buf, sizeof(buf), "%.6e\n",
                    i < data->size() ? (*data)[i] : 0.0);
      out << buf;
    }
    out << "</DataArray>\n";
  }
  for (const auto& [name, data] : int_scalars) {
    out << "<DataArray type=\"Int32\" Name=\"" << name
        << "\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < ncl; ++i)
      out << (i < data->size() ? (*data)[i] : 0) << '\n';
    out << "</DataArray>\n";
  }
  out << "</CellData>\n";

  out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
  if (!out) throw std::runtime_error("vtu: write failed: " + path);
}

}  // namespace cp
