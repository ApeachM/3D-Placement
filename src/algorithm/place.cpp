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
#include <fstream>
#include <memory>
#include "Chip.h"
#include <random>
#include <libkahypar.h>

namespace VLSI_backend {
void Chip::normalPlacement() {
  /* top die util setting manually in code level */
  vector<double> densities;
  densities.push_back(1.0);
  densities.push_back(1.0);
  setTargetDensity(densities);

  doInitialPlace();
}

enum{IGRAPH, KAHYPAR};
void Chip::partition_igraph() {
  clock_t time_start, total_time;
  int Die_1 = 0, Die_2 = 0;
  double resolution = 0.001;

  total_time = clock();

  igraph_t graph;
  igraph_vector_int_t membership;
  igraph_integer_t nb_clusters = 9999999;
  igraph_real_t quality;

  make_igraph(graph);

  igraph_vector_int_init(&membership, igraph_vcount(&graph));

  time_start = clock();
  while (nb_clusters > 20 && (resolution - 0.0001) > 1.0e-8) {
    igraph_community_leiden(&graph, nullptr, nullptr, resolution, 0.01, false, -1, &membership, &nb_clusters, &quality);
    printf("Leiden found %" IGRAPH_PRId " clusters using CPM (resolution parameter %f), quality is %.3f.\n", nb_clusters, resolution, quality);
    resolution -= 0.0001;
  }
  cout << "Leiden time: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
//  printf("Membership: ");
//  igraph_vector_int_print(&membership);

  for (int i = 0; i < instance_pointers_.size(); i++) {
    if (VECTOR(membership)[i] * 2 < nb_clusters) {
      instance_pointers_[i]->assignDie(1);
      Die_1++;
    } else {
      instance_pointers_[i]->assignDie(2);
      Die_2++;
    }
  }

  cout << "Total time: " << double(clock() - total_time) / CLOCKS_PER_SEC << "[s]" << endl;
  cout << "Dei_1: " << Die_1 << endl;
  cout << "Dei_2: " << Die_2 << endl;
  Die_2 = Die_1 = 0;
  evaluation(IGRAPH);

  igraph_vector_int_destroy(&membership);
  igraph_destroy(&graph);
}

void Chip::make_igraph(igraph_t &graph){
  igraph_integer_t num_of_cell{static_cast<int>(instance_pointers_.size())};
  igraph_integer_t num_of_net{static_cast<int>(net_pointers_.size())};

  igraph_rng_seed(igraph_rng_default(), 0);
  igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
  igraph_add_vertices(&graph, num_of_cell + num_of_net, nullptr);

  clock_t time_start = clock();

  vector<int> edges_list{};
  int net_index = static_cast<int>(instance_pointers_.size());
  for(Net *net : net_pointers_){
    for(Pin *pin : net->getConnectedPins()){
      if (pin->isInstancePin()) {
        edges_list.push_back(net_index);
        edges_list.push_back(pin->getInstance()->getId());
      }
    }
    net_index++;
  }

  if(!edges_list.empty()) {
    cout << "igraph connection: " << edges_list.size() / 2 << endl;
    igraph_vector_int_t edges;
    igraph_vector_int_init(&edges, static_cast<igraph_integer_t>(edges_list.size()));
    int x = 0;
    for (int i : edges_list) {
      VECTOR(edges)[x] = i;
      x++;
    }
    igraph_add_edges(&graph, &edges, nullptr);
    igraph_vector_int_destroy(&edges);
  }
  cout << "make igraph time: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
}

void Chip::evaluation(int type){
  int connection = 0;
  for(Net *net : net_pointers_){
    int die_1 = 0, die_2 = 0;
    for(Pin *pin : net->getConnectedPins()){
      if (pin->isInstancePin()) {
        if(pin->getInstance()->getDieId() == 1){
          die_1++;
        }
        else{
          die_2++;
        }
      }
    }
    if(die_1 && die_2){
      connection++;
    }
  }
  if(!type){
    cout << "igraph die to die connection: " << connection << endl;
  }
  else{
    cout << "KaHyPar die to die connection: " << connection << endl;
  }
}

void Chip::partition_kahypar(){
  int Die_1 = 0, Die_2 = 0;
  kahypar_context_t* context = kahypar_context_new();
//  kahypar_configure_context_from_file(context, "../submodule/kahypar/config/cut_rKaHyPar_sea20.ini");
  kahypar_configure_context_from_file(context, "../submodule/kahypar/config/cut_kKaHyPar_sea20.ini");

  clock_t time_start = clock();

  int net_index = 0;
  vector<int> hyperedge_list{};
  vector<int> hyperindices_list{};
  hyperindices_list.push_back(net_index);
  for(Net *net : net_pointers_){
    for(Pin *pin : net->getConnectedPins()){
      if (pin->isInstancePin()) {
        hyperedge_list.push_back(pin->getInstance()->getId());
        net_index++;
      }
    }
    hyperindices_list.push_back(net_index);
  }

  std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(hyperindices_list.size());
  for(int i = 0; i < hyperindices_list.size(); i++){
    hyperedge_indices[i] = hyperindices_list[i];
  }

  std::unique_ptr<kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<kahypar_hyperedge_id_t[]>(hyperedge_list.size());
  for(int i = 0; i < hyperedge_list.size(); i++){
    hyperedges[i] = hyperedge_list[i];
  }

  const double imbalance = 0.03;
  const kahypar_partition_id_t k = 2;
  kahypar_hyperedge_weight_t objective = 0;
  const kahypar_hypernode_id_t num_vertices = instance_pointers_.size();
  const kahypar_hyperedge_id_t num_hyperedges = hyperindices_list.size() - 1;
  std::vector<kahypar_partition_id_t> partition(num_vertices, -1);

  kahypar_partition(num_vertices, num_hyperedges, imbalance, k, nullptr, nullptr, hyperedge_indices.get(), hyperedges.get(), &objective, context, partition.data());

  for(int i = 0; i != num_vertices; ++i) {
    if (partition[i]) {
      instance_pointers_[i]->assignDie(2);
      Die_2++;
    }
    else {
      instance_pointers_[i]->assignDie(1);
      Die_1++;
    }
  }

  kahypar_context_free(context);

  cout << "KaHyPar_partition time: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
  cout << "Dei_1: " << Die_1 << endl;
  cout << "Dei_2: " << Die_2 << endl;
  evaluation(KAHYPAR);
}

void Chip::partition_kahypar_bin(){
  string kahypar_path = "../submodule/tool/KaHyPar";
  string graph_path = "../submodule/tool/graph.hgr";
  string config_path = "../submodule/kahypar/config/cut_kKaHyPar_sea20.ini";

  std::ofstream write(graph_path);

  const double imbalance = 0.03;

  clock_t time_start = clock();

  write << net_pointers_.size()-112 << " " << instance_pointers_.size() << "\n";

  int net_count = 0;
  std::set<size_t> cell_ids;
  for (Net *net : net_pointers_) {
    for (Pin *pin : net->getConnectedPins()) {
      if(pin->isInstancePin()){
        cell_ids.insert(pin->getInstance()->getId() + 1);
      }
    }
    if(!cell_ids.empty()){
      net_count++;
      for (auto cell_id : cell_ids) {
        write << cell_id << " ";
      }
      write << "\n";
    }
    cell_ids.clear();
  }
  write.close();

  std::stringstream command;

  command << kahypar_path << " -h " << graph_path << " -k 2" << " -e " << imbalance << " -o cut -m direct -p " << config_path << "\n";
//  command << kahypar_path << " -h " << graph_path << " -k 2" << " -e " << imbalance << " -o cut -m direct -p " << config_path << " -q true" << "\n";
  system(command.str().c_str());

  std::stringstream read_graph;
  read_graph << graph_path << ".part2" << ".epsilon" << imbalance << ".seed-1" << ".KaHyPar";

  int Die_1 = 0, Die_2 = 0;
  std::ifstream ifs(read_graph.str());
  assert(ifs.is_open());
  for (auto instance : instance_pointers_) {
    int die_id;
    ifs >> die_id;
    if (die_id) {
      instance->assignDie(2);
      Die_2++;
    }
    else {
      instance->assignDie(1);
      Die_1++;
    }
  }
  ifs.close();
  cout << "Dei_1: " << Die_1 << endl;
  cout << "Dei_2: " << Die_2 << endl;
  evaluation(KAHYPAR);
  cout << "KaHyPar_partition_binary time: " << double(clock() - time_start) / CLOCKS_PER_SEC << "[s]" << endl;
}

void Chip::placement2DieSynchronously() {

}

void Chip::do3DPlace() {
  // 1. do3DPlace the cells in the pseudo die
//  this->normalPlacement();

  // 2. partition
  this->partition_igraph();
//  this->partition_kahypar();
//  this->partition_kahypar_bin();

  // 3. placement synchronously
//  this->placement2DieSynchronously();
}
}