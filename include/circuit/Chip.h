///////////////////////////////////////////////////////////////////////////////
// Creator: Minjae Kim of CSDL, POSTECH
// Email:   kmj0824@postech.ac.kr
// GitHub:  ApeachM
//
// BSD 3-Clause License
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

#ifndef PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#define PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#include <utility>
#include <vector>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include "Parser.h"
#include "Instance.h"
#include "Net.h"
#include "Pin.h"
#include "Die.h"
#include "fft.h"
#include "HierRTLMP.h"
#include "Drawer.h"

#define REPLACE_SQRT2 1.414213562373095048801L

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;
typedef Eigen::Triplet<float> T;
using Eigen::BiCGSTAB;
using Eigen::IdentityPreconditioner;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;

namespace VLSI_backend {
using namespace odb;
// https://codingforspeed.com/using-faster-exponential-approximation/
static double fastExp(float a);

/**\brief
 * This will do 3d placement (D2D placement). \n
 *
 * \usage
 * You should call the methods following below order in normal way.\n
 * Core method is \c do3DPlace(). \n
 * <ol><li> \c parse() \n
 * <li> \c do3DPlace() \n
 * <li> \c write() \n
 *
 * \details
 * In the \c do3DPlace(), placement, partitioning, and synchronized placing will be conducted only one \c odb.\n
 * After all processing, two \c odb(database) will be generated when writing two defs.
 *
 * \author
 * Minjae Kim \n
 * GitHub: ApeachM (https://github.com/ApeachM)
 * */
class Chip {
  /* Placers */
  class InitialPlacer;
  class NesterovPlacer;
  class Partitioner;

 public:
  Chip();
  ~Chip() = default;
  /**
   * \brief
   * Core method
   * \details
   * 1. do replace in pseudo die \n
   * 2. partition into two dies
   * 3. do placement synchronously
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void do3DPlace(const string &def_name, const string &lef_name = "");

  void setBenchType(BENCH_TYPE bench_type);
  const string &getDesignName() const;

  void test();

  // etc
  void dbTutorial() const;

 protected:
  // Data initialization
  void dataBaseInit();

  /**
   * \brief
   * Read and Write functions
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void parse(const string &def_name, const string &lef_name = "");
  void write(const string &file_name);
  void parseNORMAL(const string &lef_name, const string &def_name);
  void writeNORMAL(const string &out_file_name);
  /**
   * \brief
   * This code parses and writes the input data of ICCAD 2022 contest.
   *
   * \details
   * This code highly refers to https://github.com/csdl-projects/ICCAD2022/blob/main/src/utils/Parser.cpp \n
   * The input file is from http://iccad-contest.org/Problems.html \n
   * The bench name is case1, case2, ... hidden_case1, hidden_case2, .... \n
   * In one bench, the information about tech info(lef info) and netlist info.
   *
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void parseICCAD(const string &input_file_name);
  void writeICCAD(const string &output_file_name);
  void parseICCADBenchData(const string &input_file_name);
  void odbConstructionForICCAD();

  void printDataInfo() const;

  /**\brief
   * One die (virtual die) placement before partitioning
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void normalPlacement();

  /**
   * \brief
   * Divide a cells into two circuit.
   * use Heir-RLTMP in OpenROAD. (this is in mpl2 in OpenROAD)
   * \todo
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void partition();
  /**
   *
   * */
  void partitionSimple();
  bool checkPartitionFile();
  void readPartitionFile();

  /**\brief
   * After partitioning, the
   * */
  void generateHybridBonds();

  /**
   * \brief
   * Do placement 2 Die synchronously.\n
   * This will consider the interaction between two die. \n
   * \details
   * This function is called in `do3DPlace()`.
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void placement2DieSynchronously();

  void updateHybridBondPositions();

  /**
 * \author
 * Minjae Kim \n
 * GitHub: ApeachM (https://github.com/ApeachM)
 * */
  void setTargetDensity(vector<double> densities);
  /**
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doInitialPlace();
  /**
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doNesterovPlace();

  /**\brief
   * get unit of micro
   * \details
   * the coordinate in this circuit is `return value`/1um.
   * \example
   * if the return value is 100, then
   * (20000, 30000) means coordinate (200um, 300um)
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  int getUnitOfMicro() const;

  /**
   * \brief
   * get HPWL of total circuit
   * \details
   * This function gets the HPWL by summation of all each nets in `net_pointers`
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  ulong getHPWL();
  /**
   * \brief
   * make log file for each net
   * \pre
   * `ulong getHPWL()`
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void checkHPWLForEachNet(int iteration);

  int getInstanceNumber() const;
  void setInstanceNumber(int instance_number);
  int getNetNumber() const;
  void setNetNumber(int net_number);
  dbDatabase *getDbDatabase() const;
  void setDbDatabase(dbDatabase *db_database);
  void setDesignName(const string &input_file_name);
  void drawTotalCircuit(const string &die_name = "die", bool high_resolution = false);
  void dbCapture(const string &file_name);
  void dbCaptureRead(const string &file_name);
  bool checkDbFile();
  void getAverageInstanceSize();
  void setTargetDensityManually();
  void partitionTriton();
  void setStartTime();

 protected:
  enum PHASE {
    START,
    INITIAL_PLACE,
    PARTITION,
    GENERATE_HYBRID_BOND,
    TWO_DIE_PLACE,
    END
  };
  /// this is for connection between the objects
  struct DataMapping {
    std::unordered_map<dbInst *, Instance *> inst_map;
    std::unordered_map<dbNet *, Net *> net_map;
    /// mapping_ for terminals on instance (pins on cell)
    std::unordered_map<dbITerm *, Pin *> pin_map_i;
    /// mapping_ for terminals on blocks (includes fixed pins on die)
    std::unordered_map<dbBTerm *, Pin *> pin_map_b;
  };
  DataMapping mapping_;

  PHASE phase_ = START;
  BENCH_TYPE bench_type_ = NORMAL;
  ICCAD2022BenchInformation bench_information_;
  string design_name_;
  utl::Logger logger_;

  // For pseudo die
  odb::dbDatabase *pseudo_db_database_{};

  // For top and bottom die.
  // This should be only used when parse ICCAD contest benchmark,
  // and write the two lef and def files for top and bottom die
  std::vector<odb::dbDatabase *> db_databases_{};

  Parser parser_;
  data_storage data_storage_;
  // DataMapping data_mapping_;

  std::vector<Instance *> instance_pointers_;
  std::vector<Net *> net_pointers_;
  std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
  std::vector<Pin *> pad_pointers_;
  std::vector<Die *> die_pointers_;
  std::vector<HybridBond *> hybrid_bond_pointers_;

  int num_technologies_ = 0;

  // hybrid size and spacing
  int hybrid_size_x_ = 0;
  int hybrid_size_y_ = 0;
  int hybrid_spacing_ = 0;

  int instance_number_ = 0;
  int net_number_ = 0;
  int average_instance_width_ = 0;
  int average_instance_height_ = 0;

  // first one is for top, the second one is for bottom.
  pair<int, int> max_utils_;
  // first one is for top, the second one is for bottom. This info will be copied at die.
  pair<RowInfo, RowInfo> row_infos_;

  HierRTLMPartition *hier_rtl_;
  Partitioner *partitioner_;

  string start_time_;
  ulong current_hpwl_;
};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
