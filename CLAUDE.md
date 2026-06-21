# Cprocess — Development Rules

## Architecture

Three-layer design. Never bypass a layer:

1. **C++ core** (`include/cprocess/`, `src/`) — all physics logic. Units: cm, keV, s, K.
2. **`proc::` namespace** (`include/cprocess/process.hpp`, `src/process.cpp`) — typed API that wraps the core. This is the single implementation shared by the text deck and Python. Every process step lives here.
3. **Python front-end** (`python/cprocess/`) — `Simulation` class with µm/keV/min/°C engineering units. Users call this layer only; they never call the text deck.

## Mandatory: Python wrapper for every C++ function

Whenever a new `proc::` function is added to `process.hpp`/`process.cpp`, you **must** also:

1. **Add a pybind11 binding** in `python/_cprocess.cpp` following the `proc_*` naming convention:
   ```cpp
   m.def("proc_<name>",
       [](SimState& st, <args...>) {
         std::ostringstream log;
         proc::<name>(st, <args...>, &log);
         return log.str();
       },
       py::arg("state"), <py::arg("...")...>);
   ```
2. **Add a `Simulation` method** in `python/cprocess/simulation.py` that:
   - Converts engineering units (µm → cm via `UM`, min → s via `MIN`, °C → K via `_celsius_to_k`)
   - Calls `_c.proc_<name>(self._st, ...)`
   - Appends the returned log string via `self._emit(...)`
   - Returns `self` (for method chaining) — **except** for query methods that return data

Do not implement process logic in `_cprocess.cpp` or `simulation.py`; only unit conversion and delegation belong there.

## Polygon / geometry convention

- All geometry arguments that are 2-D polygons use `std::vector<std::pair<double,double>>` in C++ (xy vertices, closed, winding-number test for inside/outside).
- In Python, the same polygons are passed as a list/sequence of `(x, y)` tuples in **micrometres**; the wrapper converts to cm with `UM = 1e-4`.
- An empty polygon list (`[]` / `{}`) means "blanket" (apply to the full surface).

## Testing

- Every new C++ `proc::` function gets a test in `tests/test_<feature>.cpp` or `tests/test_integration.cpp`.
- Every new Python method gets a test in `python/test_comprehensive.py`.
- Add new test executables to the `foreach(t ...)` loop in `CMakeLists.txt`.
- Run `cmake --build build -j$(nproc) && ctest --test-dir build -V` before committing.
- Copy `build/_cprocess*.so` to `python/cprocess/` and run `python3 python/test_comprehensive.py` before committing.

## Units summary

| Quantity    | Python API | C++ core |
|-------------|-----------|----------|
| Length      | µm        | cm       |
| Energy      | keV       | keV      |
| Time        | min       | s        |
| Temperature | °C        | K        |
| Dose        | cm⁻²      | cm⁻²     |
| Conc.       | cm⁻³      | cm⁻³     |

Conversion constants in `simulation.py`: `UM = 1e-4`, `NM = 1e-7`, `MIN = 60.0`.
