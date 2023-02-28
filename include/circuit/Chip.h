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

#define REPLACE_SQRT2 1.414213562373095048801L

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;
typedef Eigen::Triplet<float> T;
using Eigen::BiCGSTAB;
using Eigen::IdentityPreconditioner;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;

namespace VLSI_backend {
using namespace odb;
// https://codingforspeed.com/using-faster-exponential-approximation/
static float fastExp(float a);

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
  Chip() = default;
  ~Chip() = default;
  /*!
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

  /*!
   * \brief
   * Read and Write functions
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void parse(const string &lef_name, const string &def_name);
  void parse_iccad(const string &lef_name, const string &def_name);
  void write(const string &out_file_name);

  // etc
  void dbTutorial() const;

 protected:
  // Data initialization
  void init();

  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void setTargetDensity(vector<double> densities);
  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doInitialPlace();
  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doNestrovPlace();

  /**\brief
   * One die (virtual die) placement before partitioning
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void normalPlacement();

  /*!
   * \brief
   * Divide a cells into two circuit.
   * Louvain(actually, not louvain but ledien) clustering is implemented by igraph package
   * \todo
   * This code is just temporal code now. Meaningless but only simple partition is implemented.
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void partition();

  /*!\brief
   * After partitioning, the
   * */
  void generateHybridBonds();

  /*!
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

  /*!
   * \brief
   * get HPWL of total circuit
   * \details
   * This function gets the HPWL by summation of all each nets in `net_pointers`
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  ulong getHPWL();

  Parser parser_;
  data_storage data_storage_;
  data_mapping data_mapping_;

  std::vector<Instance *> instance_pointers_;
  std::vector<Net *> net_pointers_;
  std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
  std::vector<Pin *> pad_pointers_;
  std::vector<Die *> die_pointers_;

  int num_technologies_ = 0;
  int lib_cell_num_ = 0;
  int util_ = 100;

};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
