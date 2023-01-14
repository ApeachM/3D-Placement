#include <iostream>
#include "Circuit.h"

using namespace std;

int main() {
  string lefName = "Nangate45.lef";
  string defName = "medium01.def";
  string test_path_name = "../test/benchmarks/";
  string output_path_name = "../test/output/";

  VLSI_backend::Circuit circuit;
  circuit.parse(test_path_name + lefName, test_path_name + defName);
  circuit.place();
  circuit.write(output_path_name +"output_"+ defName);

}