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
#include <ctime>
#include "Chip.h"
#include <random>

namespace VLSI_backend {
void Chip::normalPlacement() {
  /* top die util setting manually in code level */
  vector<double> densities;
  densities.push_back(1.0);
  densities.push_back(1.0);
  setTargetDensity(densities);

  doInitialPlace();
}

void Chip::partition() {
  clock_t time_start, total_time;
  int Die_1 = 0, Die_2 = 0;
  double resolution = 0.11;

  int connection = 0;

  total_time = clock();

  igraph_t graph;
  igraph_vector_int_t membership;
  igraph_integer_t nb_clusters = 9999999;
  igraph_real_t quality;

  igraph_rng_seed(igraph_rng_default(), 0);
  igraph_integer_t num_of_vertices{static_cast<int>(instance_pointers_.size())};
  igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
  igraph_add_vertices(&graph, num_of_vertices, NULL);

  make_igraph(graph);
  igraph_vector_int_init(&membership, igraph_vcount(&graph));

  time_start = clock();
  while (nb_clusters > 10) {
    resolution -= 0.01;
    igraph_community_leiden(&graph, NULL, NULL, resolution, 0.01, 0, -1, &membership, &nb_clusters, &quality);
  }
  printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter %.2f), quality is %.3f.\n", nb_clusters, resolution, quality);
  cout << "Leiden Duration: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;

  printf("Membership: ");
  igraph_vector_int_print(&membership);

  for (int i = 0; i < num_of_vertices; i++) {
    Instance *cell = instance_pointers_[i];
    if (VECTOR(membership)[i] * 2 < nb_clusters) {
      cell->assignDie(1);
      Die_1++;
    } else {
      cell->assignDie(2);
      Die_2++;
    }
  }

  igraph_vector_int_destroy(&membership);
  igraph_destroy(&graph);

  cout << "Total Duration: " << double(clock() - total_time) / CLOCKS_PER_SEC << "[s]" << endl;
  cout << "Dei_1: " << Die_1 << endl;
  cout << "Dei_2: " << Die_2 << endl;

  return;
}

void Chip::make_igraph(igraph_t &graph) {
  int connection = 0;

  clock_t time_start = clock();

//  for (Net *net : net_pointers_) {
//    vector<int> cell_indices{};
//    for (Pin *pin : net->getConnectedPins()) {
//      if (pin->isInstancePin()) {
//        cell_indices.push_back(pin->getInstance()->getId());
//      }
//    }
//
//    for (int i = 0; i < cell_indices.size() - 1; i++) {
//      for (int j = i + 1; j < cell_indices.size(); j++) {
//        igraph_integer_t START(cell_indices[i]);
//        igraph_integer_t FINAL(cell_indices[j]);
//        igraph_add_edge(&graph, START, FINAL);
//        connection++;
//      }
//    }
//  }

  for (Instance *cell : instance_pointers_) {
    for (Pin *pin : cell->getPins()) {
      if(pin->isInstancePin()){
        Net* instance_net = pin->getNet();
        if(instance_net) {
          for(Pin *p : instance_net->getConnectedPins()){
            if(p->isInstancePin() && cell->getId() < p->getInstance()->getId()){
              igraph_add_edge(&graph, cell->getId(), p->getInstance()->getId());
              connection++;
            }
          }
        }
      }
    }
  }
//  cout << "Connection: " << connection << endl;
  cout << "make_igraph Duration: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
}

void Chip::placement2DieSynchronously() {

}

void Chip::do3DPlace() {
  // 1. do3DPlace the cells in the pseudo die
//  this->normalPlacement();

  // 2. partition
  this->partition();

  // 3. placement synchronously
//  this->placement2DieSynchronously();
}
}
