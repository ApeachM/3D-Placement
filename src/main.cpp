#include <iostream>
#include <igraph.h>
#include "D2DChip.h"

using namespace std;

int main() {
  string lefName = "Nangate45.lef";
  string defName = "simple01.def";
  string test_path_name = "../test/benchmarks/etc/";
  string output_path_name = "../test/output/";

  VLSI_backend::D2DChip chip;
  chip.parse(test_path_name + lefName, test_path_name + defName);
  chip.partition();
}