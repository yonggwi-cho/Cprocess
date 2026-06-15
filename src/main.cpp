#include <fstream>
#include <iostream>
#include <string>

#include "cprocess/deck.hpp"

static const char* kVersion = "cprocess 0.1.0";

int main(int argc, char** argv) {
  if (argc >= 2 && std::string(argv[1]) == "--version") {
    std::cout << kVersion << "\n";
    return 0;
  }
  if (argc < 2) {
    std::cerr << kVersion
              << "\n3D semiconductor process simulator "
                 "(unstructured FVM: implant + diffusion)\n\n"
                 "usage: cprocess <deck-file>\n"
                 "       cprocess -            (read deck from stdin)\n";
    return 2;
  }

  cp::SimState st;
  try {
    if (std::string(argv[1]) == "-") {
      cp::run_deck(std::cin, st, std::cout);
    } else {
      std::ifstream f(argv[1]);
      if (!f) {
        std::cerr << "Error: cannot open deck file: " << argv[1] << "\n";
        return 1;
      }
      cp::run_deck(f, st, std::cout);
    }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
