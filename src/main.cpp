#include <fstream>
#include <iostream>
#include <string>

#ifdef CPROCESS_MPI
#  include <mpi.h>
#endif

#include "cprocess/deck.hpp"

static const char* kVersion = "cprocess 0.1.0";

int main(int argc, char** argv) {
#ifdef CPROCESS_MPI
  // PA-5: MPI is initialized/finalized by the CLI entry point (launched via
  // mpirun), not by the library itself. cprocess_core only checks
  // MPI_Initialized and falls back to rank=0/nprocs=1 if this is not called
  // (e.g. the library is embedded in a non-MPI host application).
  MPI_Init(&argc, &argv);
#endif

  int rc = 0;
  if (argc >= 2 && std::string(argv[1]) == "--version") {
    std::cout << kVersion << "\n";
  } else if (argc < 2) {
    std::cerr << kVersion
              << "\n3D semiconductor process simulator "
                 "(unstructured FVM: implant + diffusion)\n\n"
                 "usage: cprocess <deck-file>\n"
                 "       cprocess -            (read deck from stdin)\n";
    rc = 2;
  } else {
    cp::SimState st;
    try {
      if (std::string(argv[1]) == "-") {
        cp::run_deck(std::cin, st, std::cout);
      } else {
        std::ifstream f(argv[1]);
        if (!f) {
          std::cerr << "Error: cannot open deck file: " << argv[1] << "\n";
          rc = 1;
        } else {
          cp::run_deck(f, st, std::cout);
        }
      }
    } catch (const std::exception& e) {
      std::cerr << "Error: " << e.what() << "\n";
      rc = 1;
    }
  }

#ifdef CPROCESS_MPI
  MPI_Finalize();
#endif
  return rc;
}
