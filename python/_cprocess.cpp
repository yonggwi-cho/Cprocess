#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

#include "cprocess/deck.hpp"
#include "cprocess/diffusion.hpp"
#include "cprocess/implant.hpp"
#include "cprocess/materials.hpp"
#include "cprocess/mc_implant.hpp"
#include "cprocess/mesh.hpp"
#include "cprocess/process.hpp"
#include "cprocess/vtk_writer.hpp"

namespace py = pybind11;
using namespace cp;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static py::array_t<double> vec_to_np(const std::vector<double>& v) {
  py::array_t<double> a(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

static py::array_t<double> cents_to_np(const std::vector<Vec3>& c) {
  const py::ssize_t n = static_cast<py::ssize_t>(c.size());
  py::array_t<double> a({n, py::ssize_t(3)});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < n; ++i) {
    r(i, 0) = c[i].x;
    r(i, 1) = c[i].y;
    r(i, 2) = c[i].z;
  }
  return a;
}

// Convert numpy bool/int8/uint8 array to silicon mask vector.
static std::vector<char> np_to_mask(py::array_t<std::uint8_t> arr) {
  auto r = arr.unchecked<1>();
  std::vector<char> m(static_cast<std::size_t>(r.shape(0)));
  for (py::ssize_t i = 0; i < r.shape(0); ++i)
    m[static_cast<std::size_t>(i)] = static_cast<char>(r(i) ? 1 : 0);
  return m;
}

// ---------------------------------------------------------------------------
// Module
// ---------------------------------------------------------------------------

PYBIND11_MODULE(_cprocess, m) {
  m.doc() =
      "3-D semiconductor process simulator (C++ core).\n\n"
      "Typical use::\n\n"
      "    import cprocess\n"
      "    st = cprocess.SimState()\n"
      "    log = cprocess.run_deck('''\n"
      "        mesh box xmax=0.4um ymax=0.4um zmax=0.8um nx=8 ny=8 nz=40\n"
      "        init species=P conc=1e15\n"
      "        implant species=B energy=50keV dose=1e13 method=mc ions=50000\n"
      "        diffuse time=10min temp=1000C\n"
      "    ''', st)\n"
      "    B = st.get_field('B')   # numpy array, cm^-3\n";

  // -----------------------------------------------------------------------
  // Mesh (read-only view; mesh lives inside SimState or is returned by
  // make_box_mesh; not intended to be constructed by hand from Python)
  // -----------------------------------------------------------------------
  py::class_<Mesh>(m, "Mesh")
      .def_property_readonly("n_cells",
          [](const Mesh& me) { return static_cast<int>(me.cells.size()); })
      .def_property_readonly("n_nodes",
          [](const Mesh& me) { return static_cast<int>(me.nodes.size()); })
      .def_property_readonly("cell_vol",
          [](const Mesh& me) { return vec_to_np(me.cell_vol); },
          "Per-cell volume [cm^3], shape (n_cells,)")
      .def_property_readonly("cell_cent",
          [](const Mesh& me) { return cents_to_np(me.cell_cent); },
          "Per-cell centroid [cm], shape (n_cells, 3)")
      .def_property_readonly("nodes",
          [](const Mesh& me) {
            const py::ssize_t n = static_cast<py::ssize_t>(me.nodes.size());
            py::array_t<double> a({n, py::ssize_t(3)});
            auto r = a.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < n; ++i) {
              r(i, 0) = me.nodes[i].x;
              r(i, 1) = me.nodes[i].y;
              r(i, 2) = me.nodes[i].z;
            }
            return a;
          },
          "Node coordinates [cm], shape (n_nodes, 3)")
      .def_property_readonly("cells",
          [](const Mesh& me) {
            const py::ssize_t n = static_cast<py::ssize_t>(me.cells.size());
            py::array_t<int> a({n, py::ssize_t(4)});
            auto r = a.mutable_unchecked<2>();
            for (py::ssize_t i = 0; i < n; ++i)
              for (int k = 0; k < 4; ++k) r(i, k) = me.cells[i][k];
            return a;
          },
          "Tet connectivity (node indices), shape (n_cells, 4)");

  m.def("make_box_mesh",
      [](double x0, double x1, double y0, double y1,
         double z0, double z1, int nx, int ny, int nz) {
        return make_box_mesh(x0, x1, y0, y1, z0, z1, nx, ny, nz);
      },
      py::arg("x0"), py::arg("x1"),
      py::arg("y0"), py::arg("y1"),
      py::arg("z0"), py::arg("z1"),
      py::arg("nx"), py::arg("ny"), py::arg("nz"),
      "Create a box mesh with Kuhn hex→tet subdivision.\n"
      "Coordinates are in cm; use e.g. x1=0.4e-4 for 0.4 µm.");

  // -----------------------------------------------------------------------
  // SimState
  // -----------------------------------------------------------------------
  py::class_<SimState>(m, "SimState",
      "Simulation state: mesh + per-species dopant fields.")
      .def(py::init<>())
      .def_property_readonly("has_mesh", [](const SimState& s) { return s.has_mesh; })
      .def_property_readonly("mesh",
          [](SimState& s) -> Mesh& { return s.mesh; },
          py::return_value_policy::reference_internal)
      .def("field_names",
          [](const SimState& s) {
            std::vector<std::string> names;
            for (const auto& [k, v] : s.fields) names.push_back(k);
            return names;
          },
          "Return list of dopant species names with non-empty fields.")
      .def("get_field",
          [](const SimState& s, const std::string& name) {
            auto it = s.fields.find(name);
            if (it == s.fields.end())
              throw py::key_error("no field: " + name);
            return vec_to_np(it->second);
          },
          py::arg("species"),
          "Get concentration field [cm^-3] as numpy array, shape (n_cells,).")
      .def("set_field",
          [](SimState& s, const std::string& name, py::array_t<double> arr) {
            auto r = arr.unchecked<1>();
            auto& v = s.fields[name];
            v.resize(static_cast<std::size_t>(r.shape(0)));
            for (py::ssize_t i = 0; i < r.shape(0); ++i)
              v[static_cast<std::size_t>(i)] = r(i);
          },
          py::arg("species"), py::arg("conc"),
          "Overwrite a concentration field from a numpy array.");

  // -----------------------------------------------------------------------
  // Deck runner
  // -----------------------------------------------------------------------
  m.def("run_deck",
      [](const std::string& text, SimState& st) {
        std::istringstream iss(text);
        std::ostringstream log;
        run_deck(iss, st, log);
        return log.str();
      },
      py::arg("text"), py::arg("state"),
      "Run a process deck (string). Returns the simulator log as a string.");

  m.def("run_deck_file",
      [](const std::string& path, SimState& st) {
        std::ifstream f(path);
        if (!f) throw std::runtime_error("cannot open file: " + path);
        std::ostringstream log;
        run_deck(f, st, log);
        return log.str();
      },
      py::arg("path"), py::arg("state"),
      "Run a process deck from a file. Returns the simulator log as a string.");

  // -----------------------------------------------------------------------
  // Analytic Gaussian implant
  // -----------------------------------------------------------------------
  py::class_<ImplantParams>(m, "ImplantParams",
      "Parameters for analytic (Gaussian) ion implantation.")
      .def(py::init<>())
      .def_readwrite("dose", &ImplantParams::dose, "Implant dose [cm^-2]")
      .def_readwrite("rp",   &ImplantParams::rp,   "Projected range [cm]")
      .def_readwrite("drp",  &ImplantParams::drp,  "Range straggle [cm]")
      .def_readwrite("drl",  &ImplantParams::drl,  "Lateral straggle [cm]")
      .def_readwrite("has_window", &ImplantParams::has_window)
      .def_readwrite("x1", &ImplantParams::x1)
      .def_readwrite("x2", &ImplantParams::x2)
      .def_readwrite("y1", &ImplantParams::y1)
      .def_readwrite("y2", &ImplantParams::y2)
      .def_readwrite("gamma", &ImplantParams::gamma)
      .def_readwrite("beta",  &ImplantParams::beta)
      .def_property("profile",
          [](const ImplantParams& p) {
            return p.profile == ImplantParams::Profile::pearson4
                       ? "pearson" : "gauss"; },
          [](ImplantParams& p, const std::string& s) {
            if (s == "pearson") p.profile = ImplantParams::Profile::pearson4;
            else if (s == "gauss") p.profile = ImplantParams::Profile::gauss;
            else throw std::runtime_error("profile must be gauss|pearson");
          })
      .def("set_dopant",
          [](ImplantParams& p, const std::string& name) {
            p.dopant = find_dopant(name);
            if (!p.dopant) throw std::runtime_error("unknown dopant: " + name);
          },
          py::arg("species"));

  m.def("apply_implant",
      [](const Mesh& mesh, py::array_t<std::uint8_t> mask,
         const ImplantParams& p) {
        auto mv = np_to_mask(mask);
        std::vector<double> conc(mesh.cells.size(), 0.0);
        double total = apply_implant(mesh, mv, p, conc);
        return py::make_tuple(total, vec_to_np(conc));
      },
      py::arg("mesh"), py::arg("silicon_mask"), py::arg("params"),
      "Apply Gaussian implant. Returns (total_atoms, conc_array).");

  // -----------------------------------------------------------------------
  // Monte Carlo BCA implant
  // -----------------------------------------------------------------------
  py::class_<McImplantParams>(m, "McImplantParams",
      "Parameters for Monte Carlo (BCA) ion implantation.")
      .def(py::init<>())
      .def_readwrite("dose",        &McImplantParams::dose,        "[cm^-2]")
      .def_readwrite("energy_kev",  &McImplantParams::energy_kev,  "[keV]")
      .def_readwrite("tilt_deg",    &McImplantParams::tilt_deg,    "Beam tilt [deg]")
      .def_readwrite("rotation_deg",&McImplantParams::rotation_deg,"Azimuth [deg]")
      .def_readwrite("ions",        &McImplantParams::ions,        "Simulated ion count")
      .def_readwrite("threads",     &McImplantParams::threads,     "OpenMP threads (0=auto)")
      .def_readwrite("seed",        &McImplantParams::seed,        "RNG seed")
      .def_readwrite("channeling",  &McImplantParams::channeling,  "Crystal channeling")
      .def_readwrite("has_window",  &McImplantParams::has_window)
      .def_readwrite("x1", &McImplantParams::x1)
      .def_readwrite("x2", &McImplantParams::x2)
      .def_readwrite("y1", &McImplantParams::y1)
      .def_readwrite("y2", &McImplantParams::y2)
      .def("set_dopant",
          [](McImplantParams& p, const std::string& name) {
            p.dopant = find_dopant(name);
            if (!p.dopant) throw std::runtime_error("unknown dopant: " + name);
          },
          py::arg("species"));

  py::class_<McImplantStats>(m, "McImplantStats")
      .def_readonly("deposited",     &McImplantStats::deposited)
      .def_readonly("backscattered", &McImplantStats::backscattered)
      .def_readonly("transmitted",   &McImplantStats::transmitted)
      .def_readonly("out_of_domain", &McImplantStats::out_of_domain)
      .def_readonly("in_mask",       &McImplantStats::in_mask)
      .def_readonly("unbinned",      &McImplantStats::unbinned)
      .def_readonly("weight",        &McImplantStats::weight)
      .def_readonly("rp",            &McImplantStats::rp,  "[cm]")
      .def_readonly("drp",           &McImplantStats::drp, "[cm]")
      .def("__repr__",
          [](const McImplantStats& s) {
            char buf[200];
            std::snprintf(buf, sizeof(buf),
                "McImplantStats(deposited=%lld, Rp=%.1fnm, dRp=%.1fnm)",
                s.deposited, s.rp * 1e7, s.drp * 1e7);
            return std::string(buf);
          });

  m.def("apply_mc_implant",
      [](const Mesh& mesh, py::array_t<std::uint8_t> mask,
         const McImplantParams& p) {
        auto mv = np_to_mask(mask);
        std::vector<double> conc(mesh.cells.size(), 0.0);
        std::vector<double> dmg(mesh.cells.size(), 0.0);
        auto st = apply_mc_implant(mesh, mv, p, conc, &dmg);
        return py::make_tuple(st, vec_to_np(conc), vec_to_np(dmg));
      },
      py::arg("mesh"), py::arg("silicon_mask"), py::arg("params"),
      "Run MC BCA implant.\n"
      "Returns (stats, conc_array, damage_array) where conc and damage are [cm^-3].");

  // -----------------------------------------------------------------------
  // Diffusion
  // -----------------------------------------------------------------------
  py::class_<DiffuseOpts>(m, "DiffuseOpts",
      "Options for the FVM diffusion solver.")
      .def(py::init<>())
      .def_readwrite("temp",      &DiffuseOpts::temp,      "Temperature [K]")
      .def_readwrite("time",      &DiffuseOpts::time,      "Anneal time [s]")
      .def_readwrite("dt",        &DiffuseOpts::dt,        "Time step [s] (0=auto)")
      .def_readwrite("field_enh", &DiffuseOpts::field_enh, "Electric-field drift")
      .def_readwrite("nonortho",  &DiffuseOpts::nonortho,  "Non-orthogonal correction")
      .def_readwrite("verbosity", &DiffuseOpts::verbosity)
      .def_readwrite("activation", &DiffuseOpts::activation,
          "Solid-solubility clamp on charge neutrality")
      .def_readwrite("temp_profile", &DiffuseOpts::temp_profile,
          "Piecewise-linear T profile [(t_s, T_K), ...]; empty = isothermal")
      .def_readwrite("species_parallel", &DiffuseOpts::species_parallel,
          "PA-3: parallelize per-species solve in run(): "
          "0=auto, 1=on, -1=off");

  py::class_<DirichletBC>(m, "DirichletBC",
      "Fixed-concentration (Dirichlet) boundary condition.")
      .def(py::init<>())
      .def_readwrite("species", &DirichletBC::species)
      .def_readwrite("patch",   &DirichletBC::patch,
          "Patch index: 0=xmin,1=xmax,2=ymin,3=ymax,4=zmin,5=zmax")
      .def_readwrite("conc",    &DirichletBC::conc, "[cm^-3]");

  // diffuse(state, opts, bcs=[]) → log string
  m.def("diffuse",
      [](SimState& st, const DiffuseOpts& opts,
         const std::vector<DirichletBC>& bcs) {
        // Build solve_mask: silicon cells = 1
        const int nc = static_cast<int>(st.mesh.cells.size());
        std::vector<char> smask(nc, 0);
        for (int i = 0; i < nc; ++i) {
          const int reg = st.mesh.cell_region.empty() ? 0
                                                      : st.mesh.cell_region[i];
          auto it = st.region_material.find(reg);
          if (it == st.region_material.end() || it->second == "silicon")
            smask[i] = 1;
        }
        std::ostringstream log;
        DiffusionSolver solver(st.mesh, smask, &log);
        std::vector<SpeciesField> fields;
        for (auto& [sym, conc] : st.fields) {
          const Dopant* d = find_dopant(sym);
          if (d) fields.push_back({d, &conc});
        }
        solver.run(fields, bcs, opts);
        return log.str();
      },
      py::arg("state"), py::arg("opts"),
      py::arg("bcs") = std::vector<DirichletBC>{},
      "Run the FVM diffusion solver on all species in the SimState.\n"
      "Returns the solver log as a string.");

  // -----------------------------------------------------------------------
  // VTK output
  // -----------------------------------------------------------------------
  m.def("save_vtu",
      [](const std::string& path, const SimState& st) {
        std::vector<std::pair<std::string, const std::vector<double>*>> scalars;
        for (const auto& [sym, conc] : st.fields)
          scalars.push_back({sym, &conc});
        write_vtu(path, st.mesh, scalars, {});
      },
      py::arg("path"), py::arg("state"),
      "Save the mesh and all concentration fields to a VTK .vtu file.");

  m.def("save_vtu_fields",
      [](const std::string& path, const Mesh& mesh,
         const std::vector<std::string>& names,
         const std::vector<py::array_t<double>>& arrays) {
        if (names.size() != arrays.size())
          throw std::invalid_argument("names and arrays must be the same length");
        std::vector<std::vector<double>> bufs;
        std::vector<std::pair<std::string, const std::vector<double>*>> scalars;
        for (std::size_t i = 0; i < names.size(); ++i) {
          auto r = arrays[i].unchecked<1>();
          bufs.emplace_back(r.shape(0));
          for (py::ssize_t j = 0; j < r.shape(0); ++j)
            bufs.back()[static_cast<std::size_t>(j)] = r(j);
          scalars.push_back({names[i], &bufs.back()});
        }
        write_vtu(path, mesh, scalars, {});
      },
      py::arg("path"), py::arg("mesh"),
      py::arg("field_names"), py::arg("field_arrays"),
      "Save arbitrary numpy fields to a VTK .vtu file.");

  // -----------------------------------------------------------------------
  // Native process API (proc::*) — drives the same logic as the text deck.
  // All lengths cm, energy keV, time s, temperature K. A captured log string
  // is returned where useful.
  // -----------------------------------------------------------------------
  m.def("proc_mesh_box",
      [](SimState& st, double x0, double x1, double y0, double y1,
         double z0, double z1, int nx, int ny, int nz) {
        std::ostringstream log;
        proc::mesh_box(st, x0, x1, y0, y1, z0, z1, nx, ny, nz, &log);
        return log.str();
      },
      py::arg("state"), py::arg("x0"), py::arg("x1"), py::arg("y0"),
      py::arg("y1"), py::arg("z0"), py::arg("z1"),
      py::arg("nx"), py::arg("ny"), py::arg("nz"));

  m.def("proc_mesh_gmsh",
      [](SimState& st, const std::string& file, double scale) {
        std::ostringstream log;
        proc::mesh_gmsh(st, file, scale, &log);
        return log.str();
      },
      py::arg("state"), py::arg("file"), py::arg("scale") = 1.0);

  m.def("proc_set_region",
      [](SimState& st, const std::string& material, int tag) {
        std::ostringstream log;
        proc::set_region(st, material, tag, &log);
        return log.str();
      },
      py::arg("state"), py::arg("material"), py::arg("tag") = -1);

  m.def("proc_resolve_region",
      [](const SimState& st, const std::string& v) {
        return proc::resolve_region(st, v);
      },
      py::arg("state"), py::arg("name_or_tag"));

  m.def("proc_init",
      [](SimState& st, const std::string& species, double conc, int region) {
        std::ostringstream log;
        proc::init(st, species, conc, region, &log);
        return log.str();
      },
      py::arg("state"), py::arg("species"), py::arg("conc"),
      py::arg("region") = -1);

  m.def("proc_implant_gauss",
      [](SimState& st, const std::string& species, double dose,
         double energy_kev, double rp, double drp, double drl,
         bool has_window, double x1, double x2, double y1, double y2,
         bool damage, const std::string& profile) {
        std::ostringstream log;
        const double atoms = proc::implant_gauss(st, species, dose, energy_kev,
            rp, drp, drl, has_window, x1, x2, y1, y2, damage, profile, &log);
        return py::make_tuple(atoms, log.str());
      },
      py::arg("state"), py::arg("species"), py::arg("dose"),
      py::arg("energy_kev") = 0.0, py::arg("rp") = 0.0, py::arg("drp") = 0.0,
      py::arg("drl") = 0.0, py::arg("has_window") = false,
      py::arg("x1") = 0.0, py::arg("x2") = 0.0, py::arg("y1") = 0.0,
      py::arg("y2") = 0.0, py::arg("damage") = false,
      py::arg("profile") = "gauss");

  // damage=True seeds the "I" field: with channeling=True from the MC's own
  // Kinchin-Pease damage array (x kFrenkelSurvival, capped at
  // kAmorphizationDensity); with channeling=False, falls back to the "+1"
  // model (no MC damage exists without channeling).
  m.def("proc_implant_mc",
      [](SimState& st, const std::string& species, double dose,
         double energy_kev, long long ions, double tilt_deg, double rotation_deg,
         unsigned long long seed, int threads, bool channeling,
         bool has_window, double x1, double x2, double y1, double y2,
         bool lateral_wrap, bool damage) {
        std::ostringstream log;
        const McImplantStats s = proc::implant_mc(st, species, dose, energy_kev,
            ions, tilt_deg, rotation_deg, seed, threads, channeling,
            has_window, x1, x2, y1, y2, lateral_wrap, damage, &log);
        return py::make_tuple(s, log.str());
      },
      py::arg("state"), py::arg("species"), py::arg("dose"),
      py::arg("energy_kev"), py::arg("ions") = 100000,
      py::arg("tilt_deg") = 0.0, py::arg("rotation_deg") = 0.0,
      py::arg("seed") = 1, py::arg("threads") = 0, py::arg("channeling") = false,
      py::arg("has_window") = false, py::arg("x1") = 0.0, py::arg("x2") = 0.0,
      py::arg("y1") = 0.0, py::arg("y2") = 0.0,
      py::arg("lateral_wrap") = false, py::arg("damage") = false);

  m.def("proc_photo",
      [](SimState& st, double thickness, int nz_add) {
        std::ostringstream log;
        proc::photo(st, thickness, nz_add, &log);
        return log.str();
      },
      py::arg("state"), py::arg("thickness"), py::arg("nz") = 4);

  m.def("proc_mask",
      [](SimState& st, double x1, double x2, double y1, double y2) {
        std::ostringstream log;
        proc::mask(st, x1, x2, y1, y2, &log);
        return log.str();
      },
      py::arg("state"), py::arg("x1"), py::arg("x2"), py::arg("y1"),
      py::arg("y2"));

  m.def("proc_mask_polygon",
      [](SimState& st,
         const std::vector<std::pair<double,double>>& poly) {
        std::ostringstream log;
        proc::mask_polygon(st, poly, &log);
        return log.str();
      },
      py::arg("state"), py::arg("poly"),
      "Open a polygon-shaped window in the photoresist.\n"
      "poly: list of (x, y) tuples in cm defining a closed polygon.\n"
      "Returns a log string.");

  m.def("proc_deposit",
      [](SimState& st, const std::string& material, double thickness,
         int nz_add,
         const std::vector<std::pair<double,double>>& poly) {
        std::ostringstream log;
        proc::deposit(st, material, thickness, nz_add, poly, &log);
        return log.str();
      },
      py::arg("state"), py::arg("material"), py::arg("thickness"),
      py::arg("nz_add") = 2,
      py::arg("poly") = std::vector<std::pair<double,double>>{},
      "Deposit a film on the top surface.\n"
      "thickness in cm; poly restricts deposition (empty = blanket).\n"
      "Returns a log string.");

  m.def("proc_etch",
      [](SimState& st, double depth,
         const std::vector<std::pair<double,double>>& poly,
         const std::string& material) {
        std::ostringstream log;
        proc::etch(st, depth, poly, material, &log);
        return log.str();
      },
      py::arg("state"), py::arg("depth"),
      py::arg("poly") = std::vector<std::pair<double,double>>{},
      py::arg("material") = std::string(""),
      "Etch the top surface down by depth (cm).\n"
      "poly restricts the etch region (empty = blanket, physically removes\n"
      "cells; non-empty retags cells as gas). material restricts the etch\n"
      "to a single material (empty = all non-gas materials).\n"
      "Returns a log string.");

  m.def("proc_etch_rate",
      [](SimState& st, const std::map<std::string, double>& rates,
         double time_s, bool isotropic,
         const std::vector<std::pair<double,double>>& poly) {
        std::ostringstream log;
        proc::etch_rate(st, rates, time_s, isotropic, poly, &log);
        return log.str();
      },
      py::arg("state"), py::arg("rates"), py::arg("time_s"),
      py::arg("isotropic") = true,
      py::arg("poly") = std::vector<std::pair<double,double>>{},
      "Level-set rate/time etch (P2-5): rates maps material name -> etch\n"
      "rate in cm/s (materials absent from the map are not attacked).\n"
      "isotropic=True undercuts under mask overhangs; isotropic=False\n"
      "removes material only where the local surface normal faces away\n"
      "from remaining solid (vertical/RIE approximation). poly restricts\n"
      "the etch to (x,y) columns inside the polygon (empty = blanket).\n"
      "Coexists with proc_etch (geometric depth/poly).\n"
      "Returns a log string.");

  m.def("proc_deposit_conformal",
      [](SimState& st, const std::string& material, double thickness_cm) {
        std::ostringstream log;
        proc::deposit_conformal(st, material, thickness_cm, &log);
        return log.str();
      },
      py::arg("state"), py::arg("material"), py::arg("thickness_cm"),
      "Conformal (isotropic level-set) deposit of `material`, `thickness_cm`\n"
      "thick measured normal-to-surface everywhere, including down\n"
      "sidewalls of an existing step/trench. Coexists with proc_deposit\n"
      "(purely vertical film growth).\n"
      "Returns a log string.");

  m.def("proc_oxidize",
      [](SimState& st, double time_s, double temp_k, bool wet) {
        std::ostringstream log;
        proc::oxidize(st, time_s, temp_k, wet, &log);
        return log.str();
      },
      py::arg("state"), py::arg("time_s"), py::arg("temp_k"),
      py::arg("wet") = false,
      "Blanket thermal oxidation (Deal-Grove). time in s, temp in K.\n"
      "Returns a log string.");

  m.def("proc_silicide",
      [](SimState& st, const std::string& metal, double temp_k,
         double time_s) {
        std::ostringstream log;
        proc::silicide(st, metal, temp_k, time_s, &log);
        return log.str();
      },
      py::arg("state"), py::arg("metal"), py::arg("temp_k"),
      py::arg("time_s"),
      "Blanket silicidation of a deposited metal film. K/s core units.\n"
      "Returns a log string.");

  m.def("proc_sper",
      [](SimState& st, double temp_k, double time_s) {
        std::ostringstream log;
        proc::sper(st, temp_k, time_s, &log);
        return log.str();
      },
      py::arg("state"), py::arg("temp_k"), py::arg("time_s"),
      "Solid-phase epitaxial regrowth (P3-c). Isothermal anneal at temp_k\n"
      "(K) for time_s (s); regrows MC-amorphized columns (fields['damage']\n"
      "above sper.amorph_density) bottom-up. Returns a log string.");

  m.def("proc_oxidize_2d",
      [](SimState& st, double time_s, double temp_k, bool wet) {
        std::ostringstream log;
        proc::oxidize_2d(st, time_s, temp_k, wet, &log);
        return log.str();
      },
      py::arg("state"), py::arg("time_s"), py::arg("temp_k"),
      py::arg("wet") = false,
      "2D/3D LOCOS oxidation (P2-4): bird's-beak lateral encroachment under\n"
      "a nitride mask, via a steady-state oxidant-diffusion solve each\n"
      "sub-step. Requires a nitride-masked region already present in the\n"
      "mesh; raises if none is found (use oxidize() for blanket oxidation).\n"
      "time in s, temp in K. Returns a log string.");

  m.def("proc_epitaxy",
      [](SimState& st, double thickness_cm, double temp_k, double time_s,
         const std::map<std::string, double>& doping, bool anneal) {
        std::ostringstream log;
        proc::epitaxy(st, thickness_cm, temp_k, time_s, doping, anneal, &log);
        return log.str();
      },
      py::arg("state"), py::arg("thickness_cm"), py::arg("temp_k"),
      py::arg("time_s"), py::arg("doping") = std::map<std::string, double>{},
      py::arg("anneal") = true,
      "Epitaxial Si growth with in-situ doping. cm/K/s core units.\n"
      "Returns a log string.");

  m.def("proc_mechanics",
      [](SimState& st, double temp_k, double dt_s) {
        std::ostringstream log;
        proc::mechanics(st, temp_k, dt_s, &log);
        return log.str();
      },
      py::arg("state"), py::arg("temp_k"), py::arg("dt_s"),
      "Linear-elastic FEM mechanics solve (P2-6): thermal-mismatch + "
      "intrinsic-film-stress eigenstrain loads, Dirichlet BCs (zmin fixed, "
      "lateral roller), cg_ilu0 solve, per-cell Voigt stress written to "
      "state fields sxx/syy/szz/sxy/syz/sxz (dyn/cm^2), then one dt_s of "
      "Maxwell relaxation. temp in K, dt in s. Returns a log string.");

  m.def("proc_refine",
      [](SimState& st, const std::string& species, double thresh, int passes,
         const std::string& axis) {
        std::ostringstream log;
        proc::refine(st, species, thresh, passes, axis, &log);
        return log.str();
      },
      py::arg("state"), py::arg("species"),
      py::arg("threshold") = 0.5, py::arg("passes") = 2,
      py::arg("axis") = std::string(""),
      "Adaptively split mesh edges where `species` has steep gradients "
      "(axis=\"\", default), or direction-aligned edges when axis is "
      "\"x\"/\"y\"/\"z\" (M-6, layered directional refinement).\n"
      "Returns a log string.");

  m.def("proc_strip",
      [](SimState& st) {
        std::ostringstream log;
        proc::strip(st, &log);
        return log.str();
      },
      py::arg("state"));

  m.def("proc_add_bc",
      [](SimState& st, const std::string& species, int patch, double conc) {
        std::ostringstream log;
        proc::add_bc(st, species, patch, conc, &log);
        return log.str();
      },
      py::arg("state"), py::arg("species"), py::arg("patch"), py::arg("conc"));

  m.def("proc_clear_bc",
      [](SimState& st) { proc::clear_bc(st, nullptr); }, py::arg("state"));

  m.def("proc_diffuse",
      [](SimState& st, const DiffuseOpts& opts) {
        std::ostringstream log;
        proc::diffuse(st, opts, &log);
        return log.str();
      },
      py::arg("state"), py::arg("opts"));

  m.def("proc_diffuse_ted",
      [](SimState& st, const DiffuseOpts& opts) {
        std::ostringstream log;
        proc::diffuse_ted(st, opts, &log);
        return log.str();
      },
      py::arg("state"), py::arg("opts"),
      "Transient enhanced diffusion: interstitial-coupled anneal using the "
      "'I' damage field seeded by implants with damage=True.");

  m.def("proc_save",
      [](SimState& st, const std::string& path) {
        std::ostringstream log;
        proc::save(st, path, &log);
        return log.str();
      },
      py::arg("state"), py::arg("path"));

  m.def("proc_save_state",
      [](SimState& st, const std::string& path) {
        std::ostringstream log;
        proc::save_state(st, path, &log);
        return log.str();
      },
      py::arg("state"), py::arg("path"));

  m.def("proc_load_state",
      [](SimState& st, const std::string& path) {
        std::ostringstream log;
        proc::load_state(st, path, &log);
        return log.str();
      },
      py::arg("state"), py::arg("path"));

  m.def("proc_export_device",
      [](SimState& st, const std::string& prefix) {
        std::ostringstream log;
        proc::export_device(st, prefix, &log);
        return log.str();
      },
      py::arg("state"), py::arg("path_prefix"));

  m.def("proc_profile1d",
      [](const SimState& st, const std::string& species, double x, double y) {
        return proc::profile1d(st, species, x, y);
      },
      py::arg("state"), py::arg("species"), py::arg("x"), py::arg("y"));

  m.def("proc_active_field",
      [](const SimState& st, const std::string& species, double temp_k) {
        return vec_to_np(proc::active_field(st, species, temp_k, nullptr));
      },
      py::arg("state"), py::arg("species"), py::arg("temp_k") = -1.0,
      "Per-cell electrically active concentration [cm^-3] (solid-solubility "
      "clamp); temp_k<=0 uses the last diffuse temperature.");

  m.def("proc_set_param",
      [](SimState& st, const std::string& key, double value) {
        std::ostringstream log;
        proc::set_param(st, key, value, &log);
        return log.str();
      },
      py::arg("state"), py::arg("key"), py::arg("value"),
      "Override a physical parameter (raw core units); see materials.hpp "
      "for the supported key set.");

  m.def("proc_get_param",
      [](const std::string& key, double fallback) {
        return proc::get_param(key, fallback);
      },
      py::arg("key"), py::arg("fallback") = 0.0);

  m.def("proc_list_params",
      []() { return proc::list_params(); },
      "Currently-set parameter overrides (not the full supported key set).");

  m.def("load_gds",
      [](const std::string& path, int layer) {
        return proc::load_gds(path, layer, nullptr);
      },
      py::arg("path"), py::arg("layer") = -1,
      "Read BOUNDARY polygons (cm) from a GDSII stream file on `layer` "
      "(layer<0: all layers). Free function -- does not take a SimState.");

  m.def("find_patch",
      [](const SimState& st, const std::string& name) {
        return st.mesh.find_patch(name);
      },
      py::arg("state"), py::arg("patch_name"),
      "Patch index for a named boundary (xmin/xmax/.../zmax), or -1.");

  // -----------------------------------------------------------------------
  // Convenience: unit conversions
  // -----------------------------------------------------------------------
  m.attr("um") = 1e-4;   // 1 µm in cm
  m.attr("nm") = 1e-7;   // 1 nm in cm
  m.attr("keV") = 1.0;   // already in keV
  m.attr("min") = 60.0;  // seconds
}
