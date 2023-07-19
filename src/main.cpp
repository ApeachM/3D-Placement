#include <iostream>
#include "Chip.h"

namespace flow3DMain {
void printUsage() {
  std::cout << "Usage: ./placer3D benchType <test_case>" << std::endl;

  std::cout << "\tExample 1:" << std::endl;
  std::cout << "\t\t./placer iccad22 case2" << std::endl;
  std::cout << "\t\t\tThis is for the case2 in the ICCAD 2022 contest" << std::endl;

  std::cout << "\tExample 2:" << std::endl;
  std::cout << "\t\t./placer iccad22 case2-1" << std::endl;
  std::cout << "\t\t\tThis is for the case2_hidden in the ICCAD 2022 contest" << std::endl;

  std::cout << "\tExample 3:" << std::endl;
  std::cout << "\t\t./placer ispd18 test1" << std::endl;
  std::cout << "\t\t\tThis is for the ispd18_test1" << std::endl;
}

bool mainParse(int argc, char **argv, string &test_case_name, string &test_dir,
               flow3D::BENCH_FORMAT &bench_format, flow3D::BENCH_TYPE &bench_type) {
  bool parse_fail = false;;
  string bench_type_string;
  if (argc != 3) {
    parse_fail = true;
  } else {
    // bench type determination from argument
    bench_type_string = argv[1];
    // test case determination from argument
    test_case_name = argv[2];

    ///////////////////////////////////////////////////////////////////////
    if (bench_type_string == "iccad22" || bench_type_string == "iccad23") {
      bench_format = flow3D::BENCH_FORMAT::ICCAD;
      if (bench_type_string == "iccad22")
        bench_type = flow3D::BENCH_TYPE::ICCAD22;
      if (bench_type_string == "iccad23")
        bench_type = flow3D::BENCH_TYPE::ICCAD23;
      test_dir = "../test/benchmarks/iccad/";
    } //////////////////////////////////////////////////////////////////////
    else if (bench_type_string == "ispd18") {
      bench_format = flow3D::BENCH_FORMAT::STANDARD;
      bench_type = flow3D::BENCH_TYPE::ISPD18;

      test_dir = "../test/benchmarks/standard/ispd/ispd18_test";
      if (test_case_name == "test1") {
        test_dir += "1/";
        test_case_name = "ispd18_test1";
      } else if (test_case_name == "test2") {
        test_dir += "2/";
        test_case_name = "ispd18_test2";
      } else if (test_case_name == "test3") {
        test_dir += "3/";
        test_case_name = "ispd18_test3";
      } else if (test_case_name == "test4") {
        test_dir += "4/";
        test_case_name = "ispd18_test4";
      } else if (test_case_name == "test5") {
        test_dir += "5/";
        test_case_name = "ispd18_test5";
      } else if (test_case_name == "test6") {
        test_dir += "6/";
        test_case_name = "ispd18_test6";
      } else if (test_case_name == "test7") {
        test_dir += "7/";
        test_case_name = "ispd18_test7";
      } else if (test_case_name == "test8") {
        test_dir += "8/";
        test_case_name = "ispd18_test8";
      } else if (test_case_name == "test9") {
        test_dir += "9/";
        test_case_name = "ispd18_test9";
      } else if (test_case_name == "test10") {
        test_dir += "10/";
        test_case_name = "ispd18_test10";
      } else {
        std::cout << "Please input the correct test_case." << std::endl << std::endl;
        parse_fail = true;
      }
    } //////////////////////////////////////////////////////////////////////
    else {
      std::cout << "Please input the correct bench_type." << std::endl << std::endl;
      parse_fail = true;
    }
  }

  if (parse_fail) {
    printUsage();
    return false;
  } else
    return true;
}

}

int main(int argc, char **argv) {
  string test_case_name, test_dir;
  flow3D::BENCH_FORMAT bench_format;
  flow3D::BENCH_TYPE bench_type;
  if (!flow3DMain::mainParse(argc, argv, test_case_name, test_dir, bench_format, bench_type))
    return 1;

  flow3D::Chip chip;
  chip.setDesign(test_case_name, test_dir, bench_format, bench_type);

  chip.do3DPlace();

  return 0;
}