#include <iostream>
#include "Chip.h"

using namespace std;
void printUsage() {
  cout << "Usage: ./placer is_contest <test_case>" << endl;
  cout << "\tExample 1:" << endl;
  cout << "\t\t./placer 1 2.0" << endl;
  cout << "\t\t\tThis is for the case2 in the ICCAD contest" << endl;
  cout << "\tExample 2:" << endl;
  cout << "\t\t./placer 1 2.1" << endl;
  cout << "\t\t\tThis is for the case2_hidden in the ICCAD contest" << endl;
  cout << "\tExample 3:" << endl;
  cout << "\t\t./placer 1 3.0" << endl;
  cout << "\t\t\tThis is for the case3 in the ICCAD contest" << endl;
  cout << "\tExample 4:" << endl;
  cout << "\t\t./placer 1 3.1" << endl;
  cout << "\t\t\tThis is for the case3_hidden in the ICCAD contest" << endl;
  cout << "\tExample 5:" << endl;
  cout << "\t\t./placer 0 simple01" << endl;
  cout << "\t\t\tThis is for the simple01.def in ispd18" << endl;
}

int mainParse(int argc,
              char **argv,
              string &test_case_name,
              string &test_dir,
              string &test_case_path,
              bool &is_for_contest) {
  if (argc != 3) {
    printUsage();
    return 1;
  } else {
    is_for_contest = atoi(argv[1]);
    if (is_for_contest) {
      test_dir = "../test/benchmarks/iccad2022/";
      test_case_name = argv[2];
      if (test_case_name == "2.0") {
        cout << "This is for the case2 in the contest." << endl;
        test_case_name = "case2.txt";
        test_case_path = test_dir + test_case_name;
      } else if (test_case_name == "2.1") {
        cout << "This is for the case2 in the contest." << endl;
        test_case_name = "case2_hidden.txt";
        test_case_path = test_dir + test_case_name;
      } else if (test_case_name == "3.0") {
        cout << "This is for the case3 in the contest." << endl;
        test_case_name = "case3.txt";
        test_case_path = test_dir + test_case_name;
      } else if (test_case_name == "3.1") {
        cout << "This is for the case3_hidden in the contest." << endl;
        test_case_name = "case3_hidden.txt";
        test_case_path = test_dir + test_case_name;
      } else if (test_case_name == "4.0") {
        cout << "This is for the case4 in the contest." << endl;
        test_case_name = "case4.txt";
        test_case_path = test_dir + test_case_name;
      } else if (test_case_name == "4.1") {
        cout << "This is for the case4_hidden in the contest." << endl;
        test_case_name = "case4_hidden.txt";
        test_case_path = test_dir + test_case_name;
      } else {
        cout << "Please input the correct test_case." << endl;
        cout << "Note: There's no case1." << endl;
        return 1;
      }
    } else {
      test_dir = "../test/benchmarks/etc/";
      test_case_name = argv[2];
      test_case_path = test_dir + test_case_name + ".def";
    }
  }
  return 0;
}

int main(int argc, char **argv) {
  bool is_for_contest;
  string test_case_name;
  string test_dir;
  string test_case_path;
  string output_path_name ;

  if (mainParse(argc, argv, test_case_name, test_dir, test_case_path, is_for_contest))
    return 1;
  if (is_for_contest) {
    VLSI_backend::Chip chip;
    chip.setBenchType("ICCAD");
    chip.do3DPlace(test_case_path);
  } else {
    VLSI_backend::Chip chip;
    chip.setBenchType("NORMAL");
    chip.do3DPlace(test_case_path, test_dir + "Nangate45.lef");
    // chip.getDbDatabase()->getChip()->getBlock()->saveLef((output_path_name + lefName).c_str());
  }
}
