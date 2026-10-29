#include "Comparison.h"

#include <TROOT.h>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
  try {
    gROOT->SetBatch(true);
    femtocompare::Options options;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--self-test") { femtocompare::RunSelfTests(); return 0; }
      if (arg == "--help") {
        std::cout << "compare_same_mixed [--file-a ROOT] [--file-b ROOT] [--task-a DIR] [--task-b DIR]\n"
                     "                   [--output DIR] [--label-a TEXT] [--label-b TEXT] [--self-test]\n";
        return 0;
      }
      if (i + 1 >= argc) throw std::runtime_error("Missing value for " + arg);
      const std::string value = argv[++i];
      if (arg == "--file-a") options.fileA = value;
      else if (arg == "--file-b") options.fileB = value;
      else if (arg == "--task-a") options.taskA = value;
      else if (arg == "--task-b") options.taskB = value;
      else if (arg == "--output") options.output = value;
      else if (arg == "--label-a") options.labelA = value;
      else if (arg == "--label-b") options.labelB = value;
      else throw std::runtime_error("Unknown argument: " + arg);
    }
    femtocompare::RunComparison(options);
    std::cout << "Results: " << std::filesystem::absolute(options.output) << '\n';
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
}
