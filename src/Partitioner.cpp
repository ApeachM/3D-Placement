#include "Partitioner.h"
#include "HierRTLMP.h"
namespace VLSI_backend {
// Hier-RTL-MP Partition //
void HierRTLMPartition::init() {
  int max_num_macro = hyper_parameters_.max_num_macro;
  int min_num_macro = hyper_parameters_.min_num_macro;
  int max_num_inst = hyper_parameters_.max_num_inst;
  int min_num_inst = hyper_parameters_.min_num_inst;
  float tolerance = hyper_parameters_.tolerance;
  int max_num_level = hyper_parameters_.max_num_level;
  float coarsening_ratio = hyper_parameters_.coarsening_ratio;
  int num_bundled_ios = hyper_parameters_.num_bundled_ios;
  int large_net_threshold = hyper_parameters_.large_net_threshold;
  int signature_net_threshold = hyper_parameters_.signature_net_threshold;
  float halo_width = hyper_parameters_.halo_width;
  float fence_lx = hyper_parameters_.fence_lx;
  float fence_ly = hyper_parameters_.fence_ly;
  float fence_ux = hyper_parameters_.fence_ux;
  float fence_uy = hyper_parameters_.fence_uy;
  float area_weight = hyper_parameters_.area_weight;
  float outline_weight = hyper_parameters_.outline_weight;
  float wirelength_weight = hyper_parameters_.wirelength_weight;
  float guidance_weight = hyper_parameters_.guidance_weight;
  float fence_weight = hyper_parameters_.fence_weight;
  float boundary_weight = hyper_parameters_.boundary_weight;
  float notch_weight = hyper_parameters_.notch_weight;
  float macro_blockage_weight = hyper_parameters_.macro_blockage_weight;
  float pin_access_th = hyper_parameters_.pin_access_th;
  float target_util = hyper_parameters_.target_util;
  float target_dead_space = hyper_parameters_.target_dead_space;
  float min_ar = hyper_parameters_.min_ar;
  int snap_layer = hyper_parameters_.snap_layer;
  std::string report_directory = hyper_parameters_.report_directory;

  setClusterSize(max_num_macro, min_num_macro, max_num_inst, min_num_inst);
  setClusterSizeTolerance(tolerance);
  setMaxNumLevel(max_num_level);
  setClusterSizeRatioPerLevel(coarsening_ratio);
  setNumBundledIOsPerBoundary(num_bundled_ios);
  setLargeNetThreshold(large_net_threshold);
  setSignatureNetThreshold(signature_net_threshold);
  setHaloWidth(halo_width);
  setGlobalFence(fence_lx, fence_ly, fence_ux, fence_uy);
  setAreaWeight(area_weight);
  setOutlineWeight(outline_weight);
  setWirelengthWeight(wirelength_weight);
  setGuidanceWeight(guidance_weight);
  setFenceWeight(fence_weight);
  setBoundaryWeight(boundary_weight);
  setNotchWeight(notch_weight);
  setMacroBlockageWeight(macro_blockage_weight);
  setPinAccessThreshold(pin_access_th);
  setTargetUtil(target_util);
  setTargetDeadSpace(target_dead_space);
  setMinAR(min_ar);
  setSnapLayer(snap_layer);
  setReportDirectory(report_directory.c_str());

  // constructSTA();
}
void HierRTLMPartition::partition() {
  // This function is mocking HierRTLMP::hierRTLMacroPlacer
  // but just until partition

  //
  // Get the database information
  //
  block_ = db_->getChip()->getBlock();
  dbu_ = db_->getTech()->getDbUnitsPerMicron();
  // report the default parameters
  logger_->report("area_weight_ = {}", area_weight_);
  logger_->report("outline_weight_ = {}", outline_weight_);
  logger_->report("wirelength_weight_ = {}", wirelength_weight_);
  logger_->report("guidance_weight_ = {}", guidance_weight_);
  logger_->report("fence_weight_ = {}", fence_weight_);
  logger_->report("boundary_weight_ = {}", boundary_weight_);
  logger_->report("notch_weight_ = {}", notch_weight_);
  logger_->report("macro_blockage_weight_ = {}", macro_blockage_weight_);
  logger_->report("halo_width_ = {}", halo_width_);

  //
  // Get the floorplan information
  //
  odb::Rect die_box = block_->getDieArea();
  floorplan_lx_ = mpl2::dbuToMicron(die_box.xMin(), dbu_);
  floorplan_ly_ = mpl2::dbuToMicron(die_box.yMin(), dbu_);
  floorplan_ux_ = mpl2::dbuToMicron(die_box.xMax(), dbu_);
  floorplan_uy_ = mpl2::dbuToMicron(die_box.yMax(), dbu_);

  odb::Rect core_box = block_->getCoreArea();
  int core_lx = static_cast<int>(mpl2::dbuToMicron(core_box.xMin(), dbu_));
  int core_ly = static_cast<int>(mpl2::dbuToMicron(core_box.yMin(), dbu_));
  int core_ux = static_cast<int>(mpl2::dbuToMicron(core_box.xMax(), dbu_));
  int core_uy = static_cast<int>(mpl2::dbuToMicron(core_box.yMax(), dbu_));

  logger_->report(
      "Floorplan Outline: ({}, {}) ({}, {}),  Core Outline: ({}, {}) ({}, {})",
      floorplan_lx_, floorplan_ly_, floorplan_ux_, floorplan_uy_,
      core_lx, core_ly, core_ux, core_uy);

  //
  // Compute metrics for dbModules
  // Associate all the hard macros with their HardMacro objects
  // and report the statistics
  //
  metrics_ = computeMetrics(block_->getTopModule());
  float core_area = (core_ux - core_lx) * (core_uy - core_ly);
  float util
      = (metrics_->getStdCellArea() + metrics_->getMacroArea()) / core_area;
  float core_util
      = metrics_->getStdCellArea() / (core_area - metrics_->getMacroArea());

  logger_->report(
      "Traversed logical hierarchy\n"
      "\tNumber of std cell instances : {}\n"
      "\tArea of std cell instances : {:.2f}\n"
      "\tNumber of macros : {}\n"
      "\tArea of macros : {:.2f}\n"
      "\tTotal area : {:.2f}\n"
      "\tDesign Utilization : {:.2f}\n"
      "\tCore Utilization: {:.2f}\n",
      metrics_->getNumStdCell(), metrics_->getStdCellArea(), metrics_->getNumMacro(), metrics_->getMacroArea(),
      metrics_->getStdCellArea() + metrics_->getMacroArea(), util, core_util);

  setDefaultThresholds();


  //
  // Initialize the physcial hierarchy tree
  // create root cluster
  //
  cluster_id_ = 0;
  // set the level of cluster to be 0
  root_cluster_ = new mpl2::Cluster(cluster_id_, std::string("root"), logger_);
  // set the design metrics as the metrics for the root cluster
  root_cluster_->addDbModule(block_->getTopModule());
  root_cluster_->setMetrics(*metrics_);
  cluster_map_[cluster_id_++] = root_cluster_;
  // assign cluster_id property of each instance
  for (auto inst : block_->getInsts()) {
    odb::dbIntProperty::create(inst, "cluster_id", cluster_id_);
  }

  //
  // model each bundled IO as a cluster under the root node
  // cluster_id:  0 to num_bundled_IOs x 4 are reserved for bundled IOs
  // following the sequence:  Left -> Top -> Right -> Bottom
  // starting at (floorplan_lx_, floorplan_ly_)
  // Map IOs to Pads if the design has pads
  //

  mapIOPads();
  createBundledIOs();

  // create data flow information
  logger_->report("\nCreate Data Flow");
  createDataFlow();

  // Create physical hierarchy tree in a post-order DFS manner
  logger_->report("\nPerform Clustering..");
  multiLevelCluster(root_cluster_);

  logger_->report("\nPrint Physical Hierachy Tree\n");
  printPhysicalHierarchyTree(root_cluster_, 0);

}
void Chip::Partitioner::init(const string &design_name) {
  // refer to the example of "OpenROAD/src/par/examples/timing-aware-partitioning"
  SetNetWeight(vector<float>{1});
  SetVertexWeight(vector<float>{1});
  SetPlacementWeight(vector<float>{});
  SetTimingParams(1, 1, 1, 1, 9.20000021e-10, true);
  SetFineTuneParams(1000, 10, 50, 1.5, 30, 9.99999975e-05, 4, 100, 10, 10, 100, 0.5, 25, true, 1, 4, 50, 1000);


  // Refer to TritonPart::PartitionDesign()
  logger_->report("========================================");
  logger_->report("[STATUS] Starting TritonPart Partitioner");
  logger_->report("========================================");
  logger_->info(utl::PAR, 168, "[INFO] Partitioning parameters**** ");
  dbu = db_->getTech()->getDbUnitsPerMicron();
  block_ = db_->getChip()->getBlock();
  // Parameters
  num_parts_ = 2;
  ub_factor_ = 2;
  seed_ = 0;
  vertex_dimensions_ = 1;  // for design partitioning, vertex weight is the area
  // of the instance
  hyperedge_dimensions_
      = 1;  // for design partitioning, hyperedge weight is the connectivity
  timing_aware_flag_ = false;
  if (timing_aware_flag_ == false) {
    top_n_ = 0;  // timing driven flow is disabled
  } else {
    top_n_ = 100000;  // extract the top_n critical timing paths
  }
  placement_flag_ = false;
  if (placement_flag_ == false) {
    placement_dimensions_ = 0;  // no placement information
  } else {
    placement_dimensions_ = 2;  // 2D canvas
  }
  fence_flag_ = false;
  if (fence_flag_ == false || fence_.IsValid() == false) {
    fence_.Reset();
  }
  // local parameters
  std::string fixed_file;
  std::string community_file;
  std::string group_file;
  solution_file = design_name + "_partition_info";
  logger_->info(utl::PAR, 102, "Number of partitions = {}", num_parts_);
  logger_->info(utl::PAR, 16, "UBfactor = {}", ub_factor_);
  logger_->info(utl::PAR, 17, "Seed = {}", seed_);
  logger_->info(utl::PAR, 18, "Vertex dimensions = {}", vertex_dimensions_);
  logger_->info(utl::PAR, 19, "Hyperedge dimensions = {}", hyperedge_dimensions_);
  logger_->info(utl::PAR, 20, "Placement dimensions = {}", placement_dimensions_);
  logger_->info(utl::PAR, 21, "Timing aware flag = {}", timing_aware_flag_);
  logger_->info(utl::PAR, 103, "Guardband flag = {}", guardband_flag_);
  logger_->info(utl::PAR, 23, "Global net threshold = {}", global_net_threshold_);
  logger_->info(utl::PAR, 24, "Top {} critical timing paths are extracted.", top_n_);
  logger_->info(utl::PAR, 25, "Fence aware flag = {}", fence_flag_);

  // set the random seed
  srand(seed_);  // set the random seed
  // read the netlist from OpenDB
  // for IO port and insts (std cells and macros),
  // there is an attribute for vertex_id
  logger_->report("========================================");
  logger_->report("[STATUS] Reading netlist**** ");
  // if the fence_flag_ is true, only consider the instances within the fence
  this->readNetList("", "", "");
  logger_->report("[STATUS] Finish reading netlist****");
}
void Chip::Partitioner::doPartitioning() {
  // call the multilevel partitioner to partition hypergraph_
  // but the evaluation is the original_hypergraph_
  MultiLevelPartition();
}
void Chip::Partitioner::writeSolution() {
  // Write out the solution.
  // Format 1: write the clustered netlist in verilog directly
  for (auto term : block_->getBTerms()) {
    auto vertex_id_property = odb::dbIntProperty::find(term, "vertex_id");
    const int vertex_id = vertex_id_property->getValue();
    if (vertex_id == -1) {
      continue;  // This instance is not used
    }
    const int partition_id = solution_[vertex_id];
    if (auto property = odb::dbIntProperty::find(term, "partition_id")) {
      property->setValue(partition_id);
    } else {
      odb::dbIntProperty::create(term, "partition_id", partition_id);
    }
  }

  for (auto inst : block_->getInsts()) {
    auto vertex_id_property = odb::dbIntProperty::find(inst, "vertex_id");
    const int vertex_id = vertex_id_property->getValue();
    if (vertex_id == -1) {
      continue;  // This instance is not used
    }
    const int partition_id = solution_[vertex_id];
    if (auto property = odb::dbIntProperty::find(inst, "partition_id")) {
      property->setValue(partition_id);
    } else {
      odb::dbIntProperty::create(inst, "partition_id", partition_id);
    }
  }

  // Format 2: write the explicit solution
  // each line :  instance_name  partition_id
  if (!solution_file.empty()) {
    std::string solution_file_name = solution_file;
    if (fence_flag_ == true) {
      // if the fence_flag_ is set to true, we need to update the solution file
      // to reflect the fence
      std::stringstream str_ss;
      str_ss.setf(std::ios::fixed);
      str_ss.precision(3);
      str_ss << ".lx_" << fence_.lx / dbu;
      str_ss << ".ly_" << fence_.ly / dbu;
      str_ss << ".ux_" << fence_.ux / dbu;
      str_ss << ".uy_" << fence_.uy / dbu;
      solution_file_name = solution_file_name + str_ss.str();
    }
    logger_->info(
        utl::PAR, 110, "Updated solution file name = {}", solution_file_name);
    std::ofstream file_output;
    file_output.open(solution_file_name);

    for (auto term : block_->getBTerms()) {
      if (auto property = odb::dbIntProperty::find(term, "partition_id")) {
        file_output << term->getName() << "  ";
        file_output << property->getValue() << "  ";
        file_output << std::endl;
      }
    }

    for (auto inst : block_->getInsts()) {
      if (auto property = odb::dbIntProperty::find(inst, "partition_id")) {
        file_output << inst->getName() << "  ";
        file_output << property->getValue() << "  ";
        file_output << std::endl;
      }
    }
    file_output.close();
  }

  logger_->report("===============================================");
  logger_->report("Exiting TritonPart");
}
void Chip::Partitioner::ReadNetlist() {
  // This function is only for the ICCAD contest
  // assign vertex_id property of each instance and each IO port
  // the vertex_id property will be removed after the partitioning
  vertex_weights_.clear();
  vertex_types_.clear();
  fixed_attr_.clear();
  community_attr_.clear();
  group_attr_.clear();
  placement_attr_.clear();
  // traverse all the instances
  int vertex_id = 0;
  // check if the fence constraint is specified

  for (auto term : block_->getBTerms()) {
    odb::dbIntProperty::create(term, "vertex_id", vertex_id++);
    vertex_types_.emplace_back(par::PORT);
    std::vector<float> vwts(vertex_dimensions_, 0.0);
    vertex_weights_.push_back(vwts);
    if (placement_flag_ == true) {
      odb::Rect box = term->getBBox();
      std::vector<float> loc{(box.xMin() + box.xMax()) / 2.0f,
                             (box.yMin() + box.yMax()) / 2.0f};
      placement_attr_.emplace_back(loc);
    }
  }
  for (auto inst : block_->getInsts()) {
    // -1 means that the instance is not used by the partitioner
    odb::dbIntProperty::create(inst, "vertex_id", -1);
    // there's no liberty cell in ICCAD contest benchmark

/*
      const sta::LibertyCell *liberty_cell = network_->libertyCell(inst);
      if (liberty_cell == nullptr) {
        continue;  // ignore the instance with no liberty
      }
*/
    odb::dbMaster *master = inst->getMaster();
    // check if the instance is a pad or a cover macro
    if (master->isPad() || master->isCover()) {
      continue;
    }
/*
    const float area = liberty_cell->area();
*/
    const float area = static_cast<float>(inst->getMaster()->getWidth() * inst->getMaster()->getHeight());
    std::vector<float> vwts(vertex_dimensions_, area);
    vertex_weights_.emplace_back(vwts);
    if (master->isBlock()) {
      vertex_types_.emplace_back(par::MACRO);
//      } else if (liberty_cell->hasSequentials()) {
//        vertex_types_.emplace_back(par::SEQ_STD_CELL);
    } else {
      vertex_types_.emplace_back(par::COMB_STD_CELL);
    }
    odb::dbIntProperty::find(inst, "vertex_id")->setValue(vertex_id++);
    if (placement_flag_ == true) {
      odb::dbBox *box = inst->getBBox();
      std::vector<float> loc{(box->xMin() + box->xMax()) / 2.0f,
                             (box->yMin() + box->yMax()) / 2.0f};
      placement_attr_.emplace_back(loc);
    }
  }

  num_vertices_ = vertex_id;

  // Check all the hyperedges,
  // we do not check the parallel hyperedges
  // because we need to consider timing graph
  hyperedges_.clear();
  hyperedge_weights_.clear();
  // Each net correponds to an hyperedge
  // Traverse the hyperedge and assign hyperedge_id to each net
  // the hyperedge_id property will be removed after partitioning
  int hyperedge_id = 0;
  for (auto net : block_->getNets()) {
    odb::dbIntProperty::create(net, "hyperedge_id", -1);
    // ignore all the power net
    if (net->getSigType().isSupply()) {
      continue;
    }
    // check the hyperedge
    int driver_id = -1;      // vertex id of the driver instance
    std::set<int> loads_id;  // vertex id of sink instances
    // check the connected instances
    for (odb::dbITerm *iterm : net->getITerms()) {
      odb::dbInst *inst = iterm->getInst();
      const int vertex_id
          = odb::dbIntProperty::find(inst, "vertex_id")->getValue();
      if (vertex_id == -1) {
        continue;  // the current instance is not used
      }
      if (iterm->getIoType() == odb::dbIoType::OUTPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }
    // check the connected IO pins
    for (odb::dbBTerm *bterm : net->getBTerms()) {
      const int vertex_id
          = odb::dbIntProperty::find(bterm, "vertex_id")->getValue();
      if (vertex_id == -1) {
        continue;  // the current bterm is not used
      }
      if (bterm->getIoType() == odb::dbIoType::INPUT) {
        driver_id = vertex_id;
      } else {
        loads_id.insert(vertex_id);
      }
    }
    // check the hyperedges
    std::vector<int> hyperedge;
    if (driver_id != -1 && !loads_id.empty()) {
      hyperedge.push_back(driver_id);
      for (auto &load_id : loads_id) {
        if (load_id != driver_id) {
          hyperedge.push_back(load_id);
        }
      }
    }
    // Ignore all the single-vertex hyperedge and large global netthreshold
    // if (hyperedge.size() > 1 && hyperedge.size() <= global_net_threshold_) {
    if (hyperedge.size() > 1) {
      hyperedges_.push_back(hyperedge);
      hyperedge_weights_.emplace_back(hyperedge_dimensions_, 1.0);
      odb::dbIntProperty::find(net, "hyperedge_id")->setValue(hyperedge_id++);
    }
  }  // finish hyperedge
  num_hyperedges_ = static_cast<int>(hyperedges_.size());

  // add timing features
  if (timing_aware_flag_ == true) {
    logger_->report("[STATUS] Extracting timing paths**** ");
    BuildTimingPaths();  // create timing paths
  }

  if (num_vertices_ == 0 || num_hyperedges_ == 0) {
    logger_->error(utl::PAR, 2677, "There is no vertices and hyperedges");
  }

  // build the timing graph
  // map each net to the timing arc in the timing graph
  std::vector<std::set<int>> hyperedges_arc_set;
  for (int e = 0; e < num_hyperedges_; e++) {
    const std::set<int> arc_set{e};
    hyperedges_arc_set.push_back(arc_set);
  }

  original_hypergraph_ = std::make_shared<par::Hypergraph>(vertex_dimensions_,
                                                           hyperedge_dimensions_,
                                                           placement_dimensions_,
                                                           hyperedges_,
                                                           vertex_weights_,
                                                           hyperedge_weights_,
                                                           fixed_attr_,
                                                           community_attr_,
                                                           placement_attr_,
                                                           vertex_types_,
                                                           hyperedge_slacks_,
                                                           hyperedges_arc_set,
                                                           timing_paths_,
                                                           logger_);
  // show the status of hypergraph
  logger_->info(utl::PAR, 174, "Netlist Information**");
  logger_->info(utl::PAR, 175, "Vertices = {}", original_hypergraph_->num_vertices_);
  logger_->info(
      utl::PAR, 176, "Hyperedges = {}", original_hypergraph_->num_hyperedges_);
  logger_->info(utl::PAR, 177, "Number of timing paths = {}", timing_paths_.size());
}

}