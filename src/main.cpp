#include <iostream>
#include "D2DChip.h"

using namespace std;

int main() {
//    VLSI_backend::D2DChip chip;
//    chip.parse_iccad2022("..//test/benchmarks/iccad2022/case1.txt");
//
  string lefName = "Nangate45.lef";
  string defName = "simple01.def";
  string test_path_name = "../test/benchmarks/";
  string output_path_name = "../output/";

  VLSI_backend::Circuit circuit;
  circuit.parse(test_path_name + lefName, test_path_name + defName);
  circuit.dbTutorial();
  circuit.place();
  circuit.write(output_path_name + "output_" + defName);

}