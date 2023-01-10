#include <iostream>
#include "Circuit.h"

using namespace std;

int main() {
  string lefName = "Nangate45.lef";
  string defName = "simple01.def";
  string test_path_name = "../test/benchmarks/";
  string output_path_name = "../output/";

  VLSI_backend::Circuit circuit;
  circuit.parse(test_path_name + lefName, test_path_name + defName);
  // circuit.dbTutorial();
  circuit.place();
  circuit.write(output_path_name +"output_"+ defName);

}