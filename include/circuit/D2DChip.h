///////////////////////////////////////////////////////////////////////////////
// Creator: Minjae Kim of CSDL, POSTECH
// Email:   kmj0824@postech.ac.kr
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

#ifndef INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_INCLUDE_DATASTRUCTURES_D2DCHIP_H_
#define INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_INCLUDE_DATASTRUCTURES_D2DCHIP_H_
#include <fstream>
//#include <igraph.h>
#include "Circuit.h"
namespace VLSI_backend {

/*!
 * `D2DChip` class includes
 * */
class D2DChip {
 private:
  int num_technologies_ = 0;
  /*!
   * \brief
   * A `Circuit` object Only for netlist
   * */
  Circuit netlist_;
  /*!
   * This `db_database_netlist_` is just for netlist before partitioning
   * */
  odb::dbDatabase *db_database_netlist_{};
  /*!
   * \brief
   * The vector includes VLSI_backend classes, one VLSI_backend is for each Tier(Die).
   * */
  vector<Circuit> circuits_;

  /*!
   * \brief
   * hybrid bond size
   * \details
   * first: x coordinate
   * second: y coordinate
   * */
  pair<int, int> hybrid_bond_size_{0, 0};
  /*!
   * \brief
   * the required spacing between 2 terminals and between terminal and die boundary
   * */
  int hybrid_bond_size_space_ = 0;

  /*!
   * \brief
   * number of instances(cells)
   * */
  int instance_num_ = 0;

  /*!
   * \brief
   * number of nets
   * */
  int net_num_ = 0;

  bool is_parsed_ = false;

  /*!
   * \brief
   * Divide a cells into two circuit.
   * Louvain(actually, not louvain but ledien) clustering is implemented by igraph package
   * */
  void partition();

 public:

  /*!
   * \brief
   * Parse Def file and Lef file
   * \details
   * First of all, this parsing is for only netlist.
   * When a placement process is required, then you should do(place cells) in the two `Circuit` obejct in `circuits_`.
   * */
  void parse(const string &lef_file_name, const string &def_file_name);

  /*!
   * \brief
   * Parsing the input of iccad2022 contest
   * \details
   * This function highly refers to https://github.com/csdl-projects/ICCAD2022/blob/main/src/utils/Parser.cpp#L6-L149 \n
   * However, this parsing function will parse correspond to the OpenDB structure
   *
   * \author Minjae (ApeachM)
   *
   * */
  void parse_iccad2022(const string &input_file_name);

  /*!
   * \brief
   * Core method
   * */
   void do3DPlace(){


   }

};

}

#endif //INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_INCLUDE_DATASTRUCTURES_D2DCHIP_H_
