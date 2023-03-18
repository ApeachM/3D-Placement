#include <iostream>
//#include <igraph.h>
#include "Chip.h"

using namespace std;

int main() {
  bool is_for_contest = true;
  if (is_for_contest) {
    VLSI_backend::Chip chip;
    chip.parseICCAD("../test/benchmarks/iccad2022/case1.txt");
    chip.do3DPlace();
    chip.writeICCAD("../test/output/case1.txt");
  } else {
    string lefName = "Nangate45.lef";
    string defName = "simple01.def";
    // string defName = "test.def";
    // string lefName = "test.lef";
    string test_path_name = "../test/benchmarks/etc/";
    string output_path_name = "../test/output/";

    VLSI_backend::Chip chip;
    chip.parse(test_path_name + lefName, test_path_name + defName);
    // chip.test();
    chip.do3DPlace();
    chip.write(output_path_name + defName);
    // chip.getDbDatabase()->getChip()->getBlock()->saveLef((output_path_name + lefName).c_str());
  }
//
}