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
#include <igraph.h>
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
  void do3DPlace();

  /**
   * \brief
   * Read and Write functions
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void parse(const string &lef_name, const string &def_name);
  void write(const string &out_file_name);

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
  void test();

  // etc
  void dbTutorial() const;

 protected:
  // Data initialization
  void init();

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

  int print_minciut(const igraph_t *graph, igraph_real_t value,
                    const igraph_vector_int_t *partition,
                    const igraph_vector_int_t *partition2,
                    const igraph_vector_int_t *cut,
                    const igraph_vector_t *capacity) {

    igraph_integer_t i, nc = igraph_vector_int_size(cut);
    igraph_bool_t directed = igraph_is_directed(graph);

    printf("mincut value: %g\n", (double) value);
    printf("first partition:  ");
    igraph_vector_int_print(partition);
    printf("second partition: ");
    igraph_vector_int_print(partition2);
    printf("edges in the cut: ");
    for (i = 0; i < nc; i++) {
      igraph_integer_t edge = VECTOR(*cut)[i];
      igraph_integer_t from = IGRAPH_FROM(graph, edge);
      igraph_integer_t to = IGRAPH_TO  (graph, edge);
      if (!directed && from > to) {
        igraph_integer_t tmp = from;
        from = to;
        to = tmp;
      }
      printf("%" IGRAPH_PRId "-%" IGRAPH_PRId " (%g), ", from, to, VECTOR(*capacity)[edge]);
    }
    printf("\n");

    return 0;
  }

  void partitionIGraph() {
    igraph_t graph;
    igraph_vector_int_t edges;
    igraph_vector_int_t partition, partition2, cut;
    igraph_vector_t weights;
    igraph_real_t value;

    igraph_vector_int_init(&partition, 0);
    igraph_vector_int_init(&partition2, 0);
    igraph_vector_int_init(&cut, 0);

    /* ----------------------------- */
    // graph construction //
    // due to hyper graph data structure, instance_number_ + net_number_.
    clock_t time_start = clock();

    igraph_rng_seed(igraph_rng_default(), 42);
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

    igraph_integer_t num_of_vertices{static_cast<int>(instance_number_ + net_number_)};
    igraph_add_vertices(&graph, num_of_vertices, nullptr);

    vector<int> edges_list{};
    vector<int> weights_list{};
    for (int net_index = 0; net_index < net_number_; ++net_index) {
      Net *net = net_pointers_.at(net_index);
      for (auto instance : net->getConnectedInstances()) {
        edges_list.push_back(net_number_ + net_index);
        edges_list.push_back(instance->getId());
        weights_list.push_back(static_cast<int>(instance->getArea()));
      }
    }
    igraph_vector_int_init(&edges, static_cast<igraph_integer_t>(edges_list.size()));
    igraph_vector_init(&weights, static_cast<igraph_integer_t>(edges_list.size() / 2));

    for (int i = 0; i < edges_list.size(); ++i) {
      VECTOR(edges)[i] = edges_list.at(i);
    }
    for (int j = 0; j < weights_list.size() / 2; ++j) {
      VECTOR(weights)[j] = weights_list.at(j);
    }
    igraph_add_vertices(&graph, num_of_vertices, nullptr);
    igraph_add_edges(&graph, &edges, nullptr);
    igraph_vector_int_destroy(&edges);

    cout << "igraph construction time: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
    /* ----------------------------- */
    igraph_mincut(&graph, &value, &partition, &partition2, &cut, &weights);
    print_minciut(&graph, value, &partition, &partition2, &cut, &weights);

    igraph_vector_destroy(&weights);
    igraph_vector_int_destroy(&partition);
    igraph_vector_int_destroy(&partition2);
    igraph_vector_int_destroy(&cut);
    igraph_destroy(&graph);
  }

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

  int getInstanceNumber() const;
  void setInstanceNumber(int instance_number);
  int getNetNumber() const;
  void setNetNumber(int net_number);
  dbDatabase *getDbDatabase() const;
  void setDbDatabase(dbDatabase *db_database);

  void drawDies(const string &die_name = "die", bool as_dot = true, bool draw_same_canvas = true);

 protected:
  utl::Logger logger_;
  // For pseudo die
  odb::dbDatabase *db_database_{};

  // For top and bottom die.
  // This should be only used when parse ICCAD contest benchmark,
  // and write the two lef and def files for top and bottom die
  std::vector<odb::dbDatabase *> db_databases_{};

  Parser parser_;
  data_storage data_storage_;
  // data_mapping data_mapping_;

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

  // first one is for top, the second one is for bottom.
  pair<int, int> max_utils_;
  // first one is for top, the second one is for bottom. This info will be copied at die.
  pair<RowInfo, RowInfo> row_infos_;

  HierRTLMPartition *hier_rtl_;
};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
