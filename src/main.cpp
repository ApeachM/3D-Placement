#include <iostream>
//#include <igraph.h>
#include "Chip.h"

using namespace std;

int main() {
//    VLSI_backend::Chip chip;
//    chip.parse_iccad2022("..//test/benchmarks/iccad2022/case1.txt");
//
  string lefName = "Nangate45.lef";
  string defName = "simple01.def";
  string test_path_name = "../test/benchmarks/etc/";
  string output_path_name = "../test/output/";

  VLSI_backend::Chip chip;
  chip.parse(test_path_name + lefName, test_path_name + defName);
  chip.do3DPlace();
  chip.write(output_path_name+defName);
}