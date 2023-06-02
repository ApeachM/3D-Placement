#include "Partitioner.h"
#include "HierRTLMP.h"
namespace VLSI_backend{
// Hier-RTL-MP partition //
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
}