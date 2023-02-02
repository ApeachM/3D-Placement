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
#include <igraph.h>

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
  double resolution = 0.03;

  total_time = clock();

  igraph_t graph;
  igraph_vector_int_t membership;
  igraph_integer_t nb_clusters;
  igraph_real_t quality;

  igraph_rng_seed(igraph_rng_default(), 0);
  igraph_integer_t num_of_vertices{static_cast<int>(instance_pointers_.size())};
  igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
  igraph_add_vertices(&graph, num_of_vertices, NULL);

  time_start = clock();
  for (Net *net : net_pointers_) {
    vector<int> cell_indices{};
    for (Pin *pin : net->getConnectedPins()) {
      if (pin->isInstancePin()) {
        cell_indices.push_back(pin->getInstance()->getId());
      }
    }

    for (int i = 0; i < cell_indices.size() - 1; i++) {
      for (int j = i + 1; j < cell_indices.size(); j++) {
        igraph_integer_t START(cell_indices[i]);
        igraph_integer_t FINAL(cell_indices[j]);
        igraph_add_edge(&graph, START, FINAL);
      }
    }
  }
  cout << "Net Partition Duration: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;

  igraph_vector_int_init(&membership, igraph_vcount(&graph));

  time_start = clock();
  igraph_community_leiden(&graph, NULL, NULL, resolution, 0.01, 0, -1, &membership, &nb_clusters, &quality);
  cout << "Leiden Duration: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;

  printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter %.3f), quality is %.4f.\n", nb_clusters, resolution, quality);
  printf("Membership: ");
  igraph_vector_int_print(&membership);

  for (int i = 0; i < instance_pointers_.size(); i++) {
    Instance *instance = instance_pointers_[i];
    instance->assignDie(VECTOR(membership)[i]);
  }

  while (false) {
    if (nb_clusters <= 15)
      break;

    igraph_vector_t weights;
    igraph_vector_init(&weights, nb_clusters);

    for (int i = 0; i < igraph_vector_int_size(&membership); i++) {
      VECTOR(weights)[VECTOR(membership)[i]]++;
    }

    printf("Weights: ");
    igraph_vector_print(&weights);

    igraph_vector_int_destroy(&membership);

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_add_vertices(&graph, nb_clusters, NULL);

    for (Net *net : net_pointers_) {
      vector<int> cell_indices{};
      for (Pin *pin : net->getConnectedPins()) {
        if (pin->isInstancePin()) {
          cell_indices.push_back(pin->getInstance()->getDieId());
        }
      }

      for (int i = 0; i < cell_indices.size() - 1; i++) {
        for (int j = i + 1; j < cell_indices.size(); j++) {
          igraph_integer_t START(cell_indices[i]);
          igraph_integer_t FINAL(cell_indices[j]);
          if (START != FINAL)
            igraph_add_edge(&graph, START, FINAL);
        }
      }
    }

    igraph_vector_int_init(&membership, igraph_vcount(&graph));
    igraph_community_leiden(&graph, &weights, NULL, 0.05, 0.01, 0, 10, &membership, &nb_clusters, &quality);

    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter 0.05), quality is %.4f.\n", nb_clusters, quality);
    printf("Membership: ");
    igraph_vector_int_print(&membership);

    for (Instance *instance : instance_pointers_) {
      instance->assignDie(VECTOR(membership)[instance->getDieId()]);
    }

    igraph_vector_destroy(&weights);
  }

  for (Instance *instance : instance_pointers_) {
    if(instance->getDieId() * 2 <  nb_clusters){
      instance->assignDie(1);
      Die_1++;
    }
    else {
      instance->assignDie(2);
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
