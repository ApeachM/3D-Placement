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

#include "D2DChip.h"

namespace VLSI_backend {
void D2DChip::parse_iccad2022(const string &input_file_name) {
}

void D2DChip::parse(const string &lef_file_name, const string &def_file_name) {
  netlist_.parse(lef_file_name, def_file_name);
}
void D2DChip::partition() {
  igraph_t graph;
  igraph_vector_int_t membership;
  igraph_vector_int_t degree;
  igraph_vector_t weights;
  igraph_integer_t nb_clusters;
  igraph_real_t quality;

  /* Set default seed to get reproducible results */
  igraph_rng_seed(igraph_rng_default(), 0);

  igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

  igraph_integer_t num_of_vertices{static_cast<int>(netlist_.getInstancePointers().size())};
  igraph_add_vertices(&graph, num_of_vertices, nullptr);

  vector<Net *> net_pointers = netlist_.getNetPointers();
  for (Net *net : net_pointers) {
    vector<int> cell_indices{};
    for (Pin *pin : net->getConnectedPins()) {
      if (pin->isInstancePin()) {
        cell_indices.push_back(pin->getInstance()->getId());
      }
    }
    for (auto i = 0; i < cell_indices.size(); ++i) {
      for (int j = 0; j < i; ++j) {
        igraph_integer_t index_from{cell_indices.at(i)};
        igraph_integer_t index_to{cell_indices.at(j)};
        igraph_add_edge(&graph, index_from, index_to);
      }
    }
  }

//  igraph_community_leiden(&graph, NULL, NULL, 0.05, 0.01, 1, 10, &membership, &nb_clusters, &quality);
//  igraph_community_leiden(&graph, /*const igraph_t *graph,*/
//                          NULL,   /*const igraph_vector_t *edge_weights,*/
//                          NULL,   /*const igraph_vector_t *node_weights,*/
//                          0.05,   /*const igraph_real_t resolution_parameter,*/
//                          0.01,   /*const igraph_real_t beta,*/
//                          10,     /*const igraph_bool_t start,*/
//                          /*-1,*/
//                          &membership,  /*igraph_vector_t *membership,*/
//                          &nb_clusters, /*igraph_integer_t *nb_clusters,*/
//                          &quality      /*igraph_real_t *quality*/
//                          );


//
//  igraph_vector_destroy(&weights);
//  igraph_vector_int_destroy(&degree);
//  igraph_vector_int_destroy(&membership);
//  igraph_destroy(&graph);
}

const Circuit &D2DChip::getNetlist() const {
  return netlist_;
}
void D2DChip::setNetlist(const Circuit &netlist) {
  netlist_ = netlist;
}
}

