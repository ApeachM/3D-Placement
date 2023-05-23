#include <iostream>
#include "gtest/gtest.h"
#include "Chip.h"

using namespace std;

int main() {
  bool is_for_contest = true;
  if (is_for_contest) {
    VLSI_backend::Chip chip;
    chip.parseICCAD("../test/benchmarks/iccad2022/case2.txt");
    chip.do3DPlace();
    // chip.writeICCAD("../test/output/case1.txt");
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
}

TEST(Contest, case2){
    VLSI_backend::Chip chip;
    chip.parseICCAD("../test/benchmarks/iccad2022/case2.txt");
    chip.do3DPlace();
    // chip.writeICCAD("../test/output/case1.txt");
}

TEST(Contest, case2_hidden){
  VLSI_backend::Chip chip;
  chip.parseICCAD("../test/benchmarks/iccad2022/case2_hidden.txt");
  chip.do3DPlace();
  // chip.writeICCAD("../test/output/case1.txt");
}
TEST(Contest, case3){
  VLSI_backend::Chip chip;
  chip.parseICCAD("../test/benchmarks/iccad2022/case3.txt");
  chip.do3DPlace();
  // chip.writeICCAD("../test/output/case1.txt");
}

TEST(Contest, case3_hidden){
  VLSI_backend::Chip chip;
  chip.parseICCAD("../test/benchmarks/iccad2022/case3_hidden.txt");
  chip.do3DPlace();
  // chip.writeICCAD("../test/output/case1.txt");
}

TEST(Contest, case4){
  VLSI_backend::Chip chip;
  chip.parseICCAD("../test/benchmarks/iccad2022/case4.txt");
  chip.do3DPlace();
  // chip.writeICCAD("../test/output/case1.txt");
}

TEST(Contest, case4_hidden){
  VLSI_backend::Chip chip;
  chip.parseICCAD("../test/benchmarks/iccad2022/case4_hidden.txt");
  chip.do3DPlace();
  // chip.writeICCAD("../test/output/case1.txt");
}