#include <utility>

#include "NesterovPlacer.h"

using namespace std;

namespace VLSI_backend {
Chip::NesterovPlacer::NesterovPlacer(odb::dbDatabase *db_database,
                                     std::vector<Instance *> instance_pointers,
                                     std::vector<Net *> net_pointers,
                                     std::vector<Pin *> pin_pointers,
                                     std::vector<Pin *> pad_pointers,
                                     Die *die_pointer) {
  db_database_ = db_database;
  // TODO: should check move whether method is proper or not.
  instance_pointers_ = std::move(instance_pointers);
  net_pointers_ = std::move(net_pointers);
  pin_pointers_ = std::move(pin_pointers);
  pad_pointers_ = std::move(pad_pointers);
  die_pointer_ = die_pointer;


  // hyper parameters
  if (instance_pointers_.size() < 1e5)
    initDensityPenalty = 0.01;
  else if (instance_pointers.size() < 1e3)
    initDensityPenalty = 0.1;

}
bool Chip::NesterovPlacer::initNestrovPlace(bool is_pseudo_die) {
  if (!is_base_initialized_) {
    // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L1054
    // bool Replace::initNesterovPlace()
    // void NesterovBase::init()
    is_base_initialized_ = true;

    if (instance_pointers_.empty()) {
      cout << "No placeable instance - skipping placement." << endl;
      return false;
    }
    setInstancesArea();

    // void NesterovBase::init()
    // Set a fixed seed
    srand(42);
    int dbu_per_micron = db_database_->getChip()->getBlock()->getDbUnitsPerMicron();

    if (is_pseudo_die) {
      for (Instance *instance : instance_pointers_) {
        // For any cell, add a random noise between -1 and 1 microns to each of its
        // x and y components. This is added to make it very unlikely that identical
        // cells connected in parallel do not start at the exact same position and
        // consequently shadow each other throughout the entire placement process
        int x_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
        int y_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
        instance->setCoordinate(instance->getCoordinate().first + x_offset,
                                instance->getCoordinate().second + y_offset);
      }
    }
    initFillerCells();
    initBins();

    // initialize fft structure based on bins
    fft_ = new gpl::FFT(bin_cnt_x_, bin_cnt_y_, bin_size_x_, bin_size_y_);

    for (Instance *instance : instance_pointers_) {
      instance->setDensityValueAsDefault();
    }
    for (Pin *pin : pin_pointers_) {
      pin->initDensityCoordinate();
    }
    updateDensitySize();
  }

  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L138
  // void NesterovPlace::init()

  initSLPStepsVars();

  for (int idx = 0; idx < instance_pointers_.size(); ++idx) {
    Instance *cell = instance_pointers_.at(idx);
    updateDensityCoordiLayoutInside(cell);
    cur_slp_coordinates_[idx] = prev_slp_coordinates_[idx] = cur_coordinates_[idx] = init_coordinates_[idx] =
        pair<float, float>{cell->getDensityCenterX(), cell->getDensityCenterY()};
  }
  // bin
  updateGCellDensityCenterLocation(cur_slp_coordinates_);
  prev_hpwl_ = getHpwl();
  cout << "[replace] np init: InitialHPWL: " << prev_hpwl_ << endl;
  // FFT update
  updateDensityForceBin();
  base_wire_length_coefficient_ = initWireLengthCoef / (static_cast<float>(bin_size_x_ + bin_size_y_) * 0.5);
  cout << "[replace] np init: BaseWireLengthCoef: " << base_wire_length_coefficient_ << endl;
  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());
  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());
  cout << "[replace] np init: OverflowArea: " << overflow_area_ << endl;
  cout << "[replace] np init: NesterovInstArea: " << nesterovInstsArea() << endl;
  cout << "[replace] np init: InitSumOverflow: " << sum_overflow_unscaled_ << endl;
  updateWireLengthCoef(sum_overflow_);
  cout << "[replace] np init: WireLengthCoef: " << wire_length_coefficient_x_ << endl;

  // WL update
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  // fill in cur_slp_sum_grads_, cur_slp_wire_length_grads_, cur_slp_density_grads_
  updateGradients(cur_slp_sum_grads_, cur_slp_wire_length_grads_, cur_slp_density_grads_);

  if (is_diverged_) { return false; }

  // approximately fill in prev_slp_coordinates_ to calculate lc vars
  updateInitialPrevSLPCoordi();

  // bin, FFT, wlen update with prevSLPCoordi.
  // prev_slp_coordinates_ 얘 안 바뀌였나..?
  updateGCellDensityCenterLocation(prev_slp_coordinates_);
  updateDensityForceBin();
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  // update previSumGrads_, prev_slp_wire_length_grads_, prev_slp_density_grads_
  updateGradients(prev_slp_sum_grads_, prev_slp_wire_length_grads_, prev_slp_density_grads_);

  if (is_diverged_) { return false; }

  cout << "[replace] np init: WireLengthGradSum: " << wire_length_grad_sum_ << endl;
  cout << "[replace] np init: DensityGradSum: " << density_grad_sum_ << endl;
  density_penalty_ = (wire_length_grad_sum_ / density_grad_sum_) * initDensityPenalty;
  cout << "[replace] np init: InitDensityPenalty: " << density_penalty_ << endl;
  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());
  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());
  cout << "[replace] np init: PrevSumOverflow: " << sum_overflow_unscaled_ << endl;
  step_length_ = getStepLength();
  cout << "[replace] np init: InitialStepLength: " << step_length_ << endl;

  if ((isnan(step_length_) || isinf(step_length_)) && recursion_cnt_init_slp_coef_ < maxRecursionInitSLPCoef) {
    initialPrevCoordiUpdateCoef *= 10;
    cout << "np init: steplength = 0 detected. Rerunning Nesterov::init() with initPrevSLPCoef: "
         << initialPrevCoordiUpdateCoef << endl;
    recursion_cnt_init_slp_coef_++;
    initNestrovPlace(is_pseudo_die);
  }

  if (isnan(step_length_) || isinf(step_length_)) {
    cout << "RePlAce diverged at initial iteration. Re-run with a smaller init_density_penalty value." << endl;
  }

  return true;
}
int Chip::NesterovPlacer::doNestrovPlace(int start_iter, bool only_one_iter) {
  // refer: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L482
  // int NesterovPlace::doNesterovPlace(int start_iter)

  // backTracking variable.
  // divergence detection
  float minSumOverflow = 1e30;
  float hpwlWithMinSumOverflow = 1e30;

  // dynamic adjustment of max_phi_coef
  bool is_max_phi_coef_changed = false;

  // snapshot saving detection
  bool isSnapshotSaved = false;

  // snapshot info
  vector<pair<float, float>> snapshot_coordinates;
  vector<pair<float, float>> snapshot_slp_coordinates;
  vector<pair<float, float>> snapshot_slp_sum_grads;

  float snapshot_a = 0;
  float snapshot_density_penalty = 0;
  float snapshot_step_length = 0;
  float snapshot_wl_coef_x = 0, snapshot_wl_coef_y = 0;

  bool is_diverge_tried_revert = false;

  // Core Nesterov Loop
  int iter = start_iter;
  for (; iter < max_nesterov_iter_; ++iter) {
    // cout << "[replace-test] np: InitSumOverflow: " << sum_overflow_unscaled_ << endl;

    prevA = curA;
    // here, prevA is a_(k), curA is a_(k+1)
    // See, the ePlace-MS paper's Algorithm 1
    curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5;
    // coeff is (a_k - 1) / ( a_(k+1) ) in paper.
    float coeff = (prevA - 1.0) / curA;

    // Back-Tracking loop
    int num_back_trak;
    for (num_back_trak = 0; num_back_trak < max_back_track_; ++num_back_trak) {
      // fill in nextCoordinates with given step_length_
      // here, the instance_pointers_ includes the fillers
      for (int k = 0; k < instance_pointers_.size(); ++k) {
        pair<float, float> next_coordinate(
            cur_slp_coordinates_[k].first + step_length_ * cur_slp_sum_grads_[k].first,
            cur_slp_coordinates_[k].second + step_length_ * cur_slp_sum_grads_[k].second);
        pair<float, float> next_slp_coordinate(
            next_coordinate.first + coeff * (next_coordinate.first - cur_coordinates_[k].first),
            next_coordinate.second + coeff * (next_coordinate.second - cur_coordinates_[k].second));
        Instance *current_instance = instance_pointers_.at(k);

        next_coordinates_.at(k) = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, next_coordinate.first),
            getDensityCoordiLayoutInsideY(current_instance, next_coordinate.second));
        next_slp_coordinates_[k] = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, next_slp_coordinate.first),
            getDensityCoordiLayoutInsideY(current_instance, next_slp_coordinate.second));
      }
      updateGCellDensityCenterLocation(next_slp_coordinates_);
      updateDensityForceBin();
      updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);
      updateGradients(next_slp_sum_grads_, next_slp_wire_length_grads_, next_slp_density_grads_);
      if (is_diverged_)
        break;
      if (stepLengthDivergeCheck())
        break;
    }

    // dynamic adjustment for
    // better convergence with
    // large designs
    if (!is_max_phi_coef_changed && sum_overflow_unscaled_ < 0.35f) {
      is_max_phi_coef_changed = true;
      max_phi_coef_ *= 0.99;
    }
    if (max_back_track_ == num_back_trak) {
      cout << "Backtracking limit reached so a small step will be taken" << endl;
    }
    if (is_diverged_) {
      break;
    }

    updateNextIter();
    printStateNesterov(iter);

    if (minSumOverflow > sum_overflow_unscaled_) {
      minSumOverflow = sum_overflow_unscaled_;
      hpwlWithMinSumOverflow = prev_hpwl_;
    }
    /*
     diverge detection on
     large max_phi_cof value + large design

     1) happen overflow < 20%
     2) Hpwl is growing
    */
    if (sum_overflow_unscaled_ < 0.3f && sum_overflow_unscaled_ - minSumOverflow >= 0.02f
        && hpwlWithMinSumOverflow * 1.2f < prev_hpwl_) {
      handleDiverge(snapshot_coordinates, snapshot_slp_coordinates, snapshot_slp_sum_grads, snapshot_a,
                    snapshot_density_penalty, snapshot_step_length, snapshot_wl_coef_x, snapshot_wl_coef_y,
                    is_diverge_tried_revert);
    }
    // if it reached target overflow
    if (finishCheck()) {
      iter = max_nesterov_iter_;
      break;
    }
    if (only_one_iter)
      return iter;
  }
  // in all case including diverge,
  // db should be updated.
  updateDB();
  if (is_diverged_) {
    cout << "log_->error(GPL, diverge_code_, diverge_msg_);" << endl;
  }
  return iter;

}
bool Chip::NesterovPlacer::finishCheck() const {
  if (sum_overflow_unscaled_ <= targetOverflow) {
    cout << "[NesterovSolve] Finished with Overflow: " << sum_overflow_unscaled_ << endl;
    return true;
  }
  return false;
}
void Chip::NesterovPlacer::printStateNesterov(int iter) const {
  if (iter == 0 || (iter + 1) % 10 == 0) {
    cout << "[NesterovSolve] Iter: " << iter + 1
         << "\toverflow: " << sum_overflow_unscaled_
         << "\tHPWL: " << prev_hpwl_
         << endl;
  }
}
bool Chip::NesterovPlacer::stepLengthDivergeCheck() {
  float newStepLength = getStepLength();

  if (isnan(newStepLength) || isinf(newStepLength)) {
    is_diverged_ = true;
    diverge_msg_ = "RePlAce diverged at newStepLength.";
    diverge_code_ = 305;
    return true;
  }
  if (newStepLength > step_length_ * 0.95) {
    step_length_ = newStepLength;
    return true;
  } else if (newStepLength < 0.01) {
    step_length_ = 0.01;
    return true;
  } else {
    step_length_ = newStepLength;
    return false;
  }
}
void Chip::NesterovPlacer::handleDiverge(const vector<pair<float, float>> &snapshotCoordi,
                                         const vector<pair<float, float>> &snapshotSLPCoordi,
                                         const vector<pair<float, float>> &snapshotSLPSumGrads,
                                         float snapshotA,
                                         float snapshotDensityPenalty,
                                         float snapshotStepLength,
                                         float snapshotWlCoefX,
                                         float snapshotWlCoefY,
                                         bool &isDivergeTriedRevert) {
  diverge_msg_ = "RePlAce divergence detected. ";
  diverge_msg_ += "Re-run with a smaller max_phi_cof value.";
  diverge_code_ = 307;
  is_diverged_ = true;


  // revert back the current density penality
  cur_coordinates_ = snapshotCoordi;
  cur_slp_coordinates_ = snapshotSLPCoordi;
  cur_slp_sum_grads_ = snapshotSLPSumGrads;
  curA = snapshotA;
  density_penalty_ = snapshotDensityPenalty;
  step_length_ = snapshotStepLength;
  wire_length_coefficient_x_ = snapshotWlCoefX;
  wire_length_coefficient_y_ = snapshotWlCoefY;

  updateGCellDensityCenterLocation(cur_coordinates_);
  updateDensityForceBin();
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  is_diverged_ = false;
  diverge_code_ = 0;
  diverge_msg_ = "";
  isDivergeTriedRevert = true;

  // turn off the RD forcely
  is_routability_need_ = false;
}

void Chip::NesterovPlacer::setInstancesArea() {
  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L798
  // void PlacerBase::init()
  for (Instance *instance : instance_pointers_) {
    if (instance->isInstance()) {
      if (instance->isFixed()) {
        // TODO
        //  we didn't specify the each die for nestrov.
        //  Therefore, we remain this part as blank
        // Check whether fixed instance is
        // within the core area
        // outside core area is none of RePlAce's business
        /*
        if(isCoreAreaOverlap){
          // ... //
        }
        */
      } else {
        uint inst_area = instance->getArea();
        place_instances_area_ += inst_area;
        // macro cells should be
        // macro_instances_area_
        uint row_height = (*db_database_->getChip()->getBlock()->getRows().begin())->getSite()->getHeight();
        if (instance->getHeight() > (row_height * 6)) {
          macro_instances_area_ += inst_area;
        }
          // smaller or equal height cells should be
          // stdInstArea_
        else {
          std_instances_area_ += inst_area;
        }
      }
    } else if (instance->isDummy()) {
      dummyInsts_.push_back(instance);
      nonPlaceInsts_.push_back(instance);
      non_place_instances_area_ += instance->getArea();
    }
  }
}
void Chip::NesterovPlacer::initFillerCells() {
  // extract average dx/dy in range (10%, 90%)
  vector<int> width_storage;
  vector<int> height_storage;
  width_storage.reserve(instance_pointers_.size());
  height_storage.reserve(instance_pointers_.size());
  for (Instance *instance : instance_pointers_) {
    width_storage.push_back(static_cast<int>(instance->getWidth()));
    height_storage.push_back(static_cast<int>(instance->getHeight()));
  }

  // sort
  std::sort(width_storage.begin(), width_storage.end());
  std::sort(height_storage.begin(), height_storage.end());

  // average from (10 - 90%)
  int64_t width_sum = 0, height_sum = 0;

  int min_idx = static_cast<int>(static_cast<float>(width_storage.size()) * 0.05);
  int max_idx = static_cast<int>(static_cast<float>(width_storage.size()) * 0.95);

  // when #instances are too small,
  // extracts average values in whole ranges.
  if (min_idx == max_idx) {
    min_idx = 0;
    max_idx = static_cast<int>(width_storage.size());
  }

  for (int i = min_idx; i < max_idx; i++) {
    width_sum += width_storage[i];
    height_sum += height_storage[i];
  }

  // the avg width and avg height will be used as filler cells
  // width and height
  filler_width_ = static_cast<int>(width_sum / (max_idx - min_idx));
  filler_height_ = static_cast<int>(height_sum / (max_idx - min_idx));

  int64_t core_area = die_pointer_->getArea();

  // nonPlaceInstsArea should not have density downscaling!!!
  white_space_area_ = core_area - non_place_instances_area_;

  if (use_uniform_target_density_) {
    target_density_ =
        static_cast<float>(std_instances_area_) / static_cast<float>(white_space_area_ - macro_instances_area_)
            + 0.01;
  } else {
    target_density_ = die_pointer_->getDensity();
  }

  // density screening
  movable_area_ = white_space_area_ * target_density_;
  total_filler_area_ = movable_area_ - nesterovInstsArea();
  uniform_target_density_ = static_cast<float>(nesterovInstsArea()) / static_cast<float>(white_space_area_);
  if (total_filler_area_ < 0) {
    uniform_target_density_ = ceilf(uniform_target_density_ * 100) / 100;
    cout << "Use a higher -density or" << endl <<
         "re-floorplan with a larger core area." << endl <<
         "Given target density: " << target_density_ << endl <<
         "Suggested target density: " << uniform_target_density_ << endl;
  }
  int fillerCnt = static_cast<int>(
      total_filler_area_ / static_cast<int64_t>(filler_width_ * filler_height_));

  mt19937 genX;
  mt19937 genY;
  if (seed_fix) {
    genX.seed(1111); // fix seed
    genY.seed(2222); // fix seed
  } else {
    genX.seed(random_device{}());
    genY.seed(random_device{}());
  }
  uniform_int_distribution<int> disX(0, (int) die_pointer_->getWidth());
  uniform_int_distribution<int> disY(0, (int) die_pointer_->getHeight());

  // make and store the fillers
  fillers_.reserve(fillerCnt);
  for (int i = 0; i < fillerCnt; ++i) {
    Instance filler;
    filler.setWidth(filler_width_);
    filler.setHeight(filler_height_);
    filler.setCoordinate(disX(genX), disY(genY));
    fillers_.push_back(filler);
  }
  // setting the pointers of fillers in real data storage
  for (auto &filler : fillers_) {
    instance_pointers_.push_back(&filler);
  }

}
void Chip::NesterovPlacer::initBins() {
  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L723
  // void BinGrid::initBins()
  int ux_ = die_pointer_->getUpperRightX();
  int lx_ = die_pointer_->getLowerLeftX();
  int uy_ = die_pointer_->getUpperRightY();
  int ly_ = die_pointer_->getLowerLeftY();
  int64_t totalBinArea = static_cast<int64_t>(ux_ - lx_) * static_cast<int64_t>(uy_ - ly_);
  int64_t averagePlaceInstArea = place_instances_area_ / (instance_pointers_.size() - fillers_.size());

  int64_t idealBinArea = std::round(static_cast<float>(averagePlaceInstArea) / target_density_);
  int idealBinCnt = totalBinArea / idealBinArea;
  if (idealBinCnt < 4)   // the smallest we allow is 2x2 bins
    idealBinCnt = 4;
/*
      log_->info(GPL, 23, "TargetDensity: {:.2f}", targetDensity_);
      log_->info(GPL, 24, "AveragePlaceInstArea: {}", averagePlaceInstArea);
      log_->info(GPL, 25, "IdealBinArea: {}", idealBinArea);
      log_->info(GPL, 26, "IdealBinCnt: {}", idealBinCnt);
      log_->info(GPL, 27, "TotalBinArea: {}", totalBinArea);
*/
  int foundBinCnt = 2;
  // find binCnt: 2, 4, 8, 16, 32, 64, ...
  // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
  for (foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2) {
    if (foundBinCnt * foundBinCnt <= idealBinCnt
        && 4 * foundBinCnt * foundBinCnt > idealBinCnt) {
      break;
    }
  }
  bin_cnt_x_ = foundBinCnt;
  bin_cnt_y_ = foundBinCnt;
  // log_->info(GPL, 28, "BinCnt: {} {}", bin_cnt_x_, bin_cnt_y_);
  bin_size_x_ = ceil(static_cast<float>((ux_ - lx_)) / bin_cnt_x_);
  bin_size_y_ = ceil(static_cast<float>((uy_ - ly_)) / bin_cnt_y_);
  // log_->info(GPL, 29, "BinSize: {} {}", bin_size_x_, bin_size_y_);
  binStor_.resize(bin_cnt_x_ * bin_cnt_y_);
  bins_.reserve(bin_cnt_x_ * bin_cnt_y_);
  int x = lx_, y = ly_;
  int idxX = 0, idxY = 0;
  for (auto &bin : binStor_) {
    int sizeX = (x + bin_size_x_ > ux_) ? ux_ - x : bin_size_x_;
    int sizeY = (y + bin_size_y_ > uy_) ? uy_ - y : bin_size_y_;
    bin = Bin(idxX, idxY, x, y, x + sizeX, y + sizeY, target_density_);

    // move x, y coordinates.
    x += bin_size_x_;
    idxX += 1;
    if (x >= ux_) {
      y += bin_size_y_;
      x = lx_;
      idxY++;
      idxX = 0;
    }
    bins_.push_back(&bin);
  }


  // only initialized once
  updateBinsNonPlaceArea();
}
void Chip::NesterovPlacer::updateBinsNonPlaceArea() {
  for (auto &bin : bins_) {
    bin->setNonPlaceArea(0);
    bin->setNonPlaceAreaUnscaled(0);
  }

  for (auto &inst : nonPlaceInsts_) {
    std::pair<int, int> pairX = getMinMaxIdxX(inst);
    std::pair<int, int> pairY = getMinMaxIdxY(inst);
    for (int i = pairX.first; i < pairX.second; i++) {
      for (int j = pairY.first; j < pairY.second; j++) {
        Bin *bin = bins_[j * bin_cnt_x_ + i];

        // Note that nonPlaceArea should have scale-down with
        // target density.
        // See MS-replace paper
        //
        bin->addNonPlaceArea(getOverlapArea(bin, inst, db_database_->getChip()->getBlock()->getDbUnitsPerMicron())
                                 * bin->targetDensity());
        bin->addNonPlaceAreaUnscaled(getOverlapAreaUnscaled(bin, inst) * bin->targetDensity());
      }
    }
  }
}
float Chip::NesterovPlacer::getDensityCoordiLayoutInsideX(Instance *instance, float cx) {
  float adjVal = cx;  // adjusted value
  // will change base on each assigned binGrids.
  if (cx - instance->getDensityDeltaX() / 2 < die_pointer_->getLowerLeftX()) {
    adjVal = die_pointer_->getLowerLeftX() + instance->getDensityDeltaX() / 2;
  }
  if (cx + instance->getDensityDeltaX() / 2 > die_pointer_->getUpperRightX()) {
    adjVal = die_pointer_->getUpperRightX() - instance->getDensityDeltaX() / 2;
  }
  return adjVal;
}
float Chip::NesterovPlacer::getDensityCoordiLayoutInsideY(Instance *instance, float cy) {
  float adjVal = cy;  // adjusted value
  // will change base on each assigned binGrids.
  if (cy - instance->getDensityDeltaY() / 2 < die_pointer_->getLowerLeftY()) {
    adjVal = die_pointer_->getLowerLeftY() + instance->getDensityDeltaY() / 2;
  }
  if (cy + instance->getDensityDeltaY() / 2 > die_pointer_->getUpperRightY()) {
    adjVal = die_pointer_->getUpperRightY() - instance->getDensityDeltaY() / 2;
  }
  return adjVal;
}
void Chip::NesterovPlacer::updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates) {
  for (int idx = 0; idx < coordinates.size(); ++idx) {
    pair<float, float> coordinate = coordinates.at(idx);
    Instance *instance = instance_pointers_.at(idx);
    instance->setDensityCenterLocation(coordinate.first, coordinate.second);
  }
  updateBinsCellDensityArea(instance_pointers_);
}
std::pair<int, int> Chip::NesterovPlacer::getMinMaxIdxX(Instance *inst) const {
  int lowerIdx = (inst->ly() - die_pointer_->getLowerLeftY()) / bin_size_y_;
  int upperIdx = (fastModulo((inst->uy() - lx()), bin_size_y_) == 0)
                 ? (inst->uy() - ly()) / bin_size_y_
                 : (inst->uy() - ly()) / bin_size_y_ + 1;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, bin_cnt_y_));
}
std::pair<int, int> Chip::NesterovPlacer::getMinMaxIdxY(Instance *inst) const {
  int lowerIdx = (inst->ly() - ly()) / bin_size_y_;
  int upperIdx = (fastModulo((inst->uy() - ly()), bin_size_y_) == 0)
                 ? (inst->uy() - ly()) / bin_size_y_
                 : (inst->uy() - ly()) / bin_size_y_ + 1;

  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, bin_cnt_y_));
}
int Chip::NesterovPlacer::fastModulo(const int input, const int ceil) { return input >= ceil ? input % ceil : input; }
int64_t Chip::NesterovPlacer::getOverlapArea(const Chip::NesterovPlacer::Bin *bin, Instance *inst, int dbu_per_micron) {
  int rectLx = max(bin->lx(), inst->lx()), rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), rectUy = min(bin->uy(), inst->uy());

  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  }

  if (inst->isMacro()) {
    float meanX = (inst->cx() - inst->lx()) / (float) dbu_per_micron;
    float meanY = (inst->cy() - inst->ly()) / (float) dbu_per_micron;

    // For the bivariate normal distribution, we are using
    // the shifted means of X and Y.
    // Sigma is used as the mean/4 for both dimensions
    const biNormalParameters i
        = {meanX,
           meanY,
           meanX / 4,
           meanY / 4,
           (rectLx - inst->lx()) / (float) dbu_per_micron,
           (rectLy - inst->ly()) / (float) dbu_per_micron,
           (rectUx - inst->lx()) / (float) dbu_per_micron,
           (rectUy - inst->ly()) / (float) dbu_per_micron};

    const float original = static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
    const float scaled = calculateBiVariateNormalCDF(i)
        * static_cast<float>(inst->ux() - inst->lx())
        * static_cast<float>(inst->uy() - inst->ly());

    // For heavily dense regions towards the center of the macro,
    // we are using an upper limit of 1.15*(overlap) between the macro
    // and the bin.
    if (scaled >= original) {
      return min<float>(scaled, original * 1.15);
    }
      // If the scaled value is smaller than the actual overlap
      // then use the original overlap value instead.
      // This is implemented to prevent cells from being placed
      // at the outer sides of the macro.
    else {
      return original;
    }
  } else {
    return static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
  }
}
float Chip::NesterovPlacer::calculateBiVariateNormalCDF(Chip::NesterovPlacer::biNormalParameters i) {
  const float x1 = (i.meanX - i.lx) / (sqrt(2) * i.sigmaX);
  const float x2 = (i.meanX - i.ux) / (sqrt(2) * i.sigmaX);

  const float y1 = (i.meanY - i.ly) / (sqrt(2) * i.sigmaY);
  const float y2 = (i.meanY - i.uy) / (sqrt(2) * i.sigmaY);

  return 0.25
      * (erf(x1) * erf(y1) + erf(x2) * erf(y2) - erf(x1) * erf(y2)
          - erf(x2) * erf(y1));
}
int64_t Chip::NesterovPlacer::getOverlapAreaUnscaled(const Chip::NesterovPlacer::Bin *bin, Instance *inst) {
  int rectLx = max(bin->lx(), inst->lx()), rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), rectUy = min(bin->uy(), inst->uy());

  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  } else {
    return static_cast<int64_t>(rectUx - rectLx)
        * static_cast<int64_t>(rectUy - rectLy);
  }
}
void Chip::NesterovPlacer::updateDensitySize() {
  for (Instance *instance : instance_pointers_) {
    float scale_x = 0, scale_y = 0;
    float density_width = 0, density_height = 0;
    if (instance->dx() < REPLACE_SQRT2 * bin_size_x_) {
      scale_x = static_cast<float>(instance->dx()) / static_cast<float>(REPLACE_SQRT2 * bin_size_x_);
      density_width = REPLACE_SQRT2 * static_cast<float>(bin_size_x_);
    } else {
      scale_x = 1.0;
      density_width = instance->dx();
    }

    if (instance->dy() < REPLACE_SQRT2 * bin_size_y_) {
      scale_y = static_cast<float>(instance->dy())
          / static_cast<float>(REPLACE_SQRT2 * bin_size_y_);
      density_height = REPLACE_SQRT2 * static_cast<float>(bin_size_y_);
    } else {
      scale_y = 1.0;
      density_height = instance->dy();
    }

    instance->setDensitySize(density_width, density_height);
    instance->setDensityScale(scale_x * scale_y);
  }
}
void Chip::NesterovPlacer::updateDensityForceBin() {
  // copy density to utilize FFT
  for (auto &bin : bins_) {
    fft_->updateDensity(bin->x(), bin->y(), bin->density());
  }

  // do FFT
  fft_->doFFT();

  // update electroPhi and electroForce
  // update sum_phi_ for nesterov loop
  sum_phi_ = 0;
  for (auto &bin : bins_) {
    auto eForcePair = fft_->getElectroForce(bin->x(), bin->y());
    bin->setElectroForce(eForcePair.first, eForcePair.second);

    float electroPhi = fft_->getElectroPhi(bin->x(), bin->y());
    bin->setElectroPhi(electroPhi);

    sum_phi_ += electroPhi * static_cast<float>(bin->nonPlaceArea() + bin->instPlacedArea() + bin->fillerArea());
  }

}
void Chip::NesterovPlacer::updateWireLengthForceWA(double wlCoeffX, double wlCoeffY) {
  // clear all WA variables.
  for (Net *gNet : net_pointers_) {
    gNet->clearWaVars();
  }
  for (auto &gPin : pin_pointers_) {
    gPin->clearWaVars();
  }

  for (Net *&gNet : net_pointers_) {
    gNet->updateBox(die_pointer_->getDieId());
    vector<Pin *> pin_set;

    if (!gNet->isIntersected())
      pin_set = gNet->getConnectedPins();
    else {
      for (Pin *pin : gNet->getConnectedPins()) {
        if (pin->isInstancePin()) {
          if (pin->getInstance()->getDieId() == die_pointer_->getDieId())
            pin_set.push_back(pin);
        } else if (pin->isBlockPin()) {
          // TODO: more accurate method is needed.
          pin_set.push_back(pin);
        } else if (pin->isHybridBondPin()) {
          pin_set.push_back(pin);
        }
      }
    }

    for (Pin *pin : pin_set) {
      // The WA terms are shift invariant:
      //
      //   Sum(x_i * exp(x_i))    Sum(x_i * exp(x_i - C))
      //   -----------------    = -----------------
      //   Sum(exp(x_i))          Sum(exp(x_i - C))
      //
      // So we shift to keep the exponential from overflowing
      double expMinX = static_cast<double>(gNet->lx() - pin->cx()) * wlCoeffX;
      double expMaxX = static_cast<double>(pin->cx() - gNet->ux()) * wlCoeffX;
      double expMinY = static_cast<double>(gNet->ly() - pin->cy()) * wlCoeffY;
      double expMaxY = static_cast<double>(pin->cy() - gNet->uy()) * wlCoeffY;

      // min x
      if (expMinX > min_wire_length_force_bar_) {
        pin->setMinExpSumX(fastExp(expMinX));
        gNet->addWaExpMinSumX(pin->minExpSumX());
        gNet->addWaXExpMinSumX(pin->cx() * pin->minExpSumX());
      }

      // max x
      if (expMaxX > min_wire_length_force_bar_) {
        pin->setMaxExpSumX(fastExp(expMaxX));
        gNet->addWaExpMaxSumX(pin->maxExpSumX());
        gNet->addWaXExpMaxSumX(pin->cx() * pin->maxExpSumX());
      }

      // min y
      if (expMinY > min_wire_length_force_bar_) {
        pin->setMinExpSumY(fastExp(expMinY));
        gNet->addWaExpMinSumY(pin->minExpSumY());
        gNet->addWaYExpMinSumY(pin->cy() * pin->minExpSumY());
      }

      // max y
      if (expMaxY > min_wire_length_force_bar_) {
        pin->setMaxExpSumY(fastExp(expMaxY));
        gNet->addWaExpMaxSumY(pin->maxExpSumY());
        gNet->addWaYExpMaxSumY(pin->cy() * pin->maxExpSumY());
      }
    }
  }

}
float Chip::NesterovPlacer::nesterovInstsArea() const {
  return std_instances_area_ + static_cast<int64_t>(round(macro_instances_area_ * target_density_));
}
void Chip::NesterovPlacer::updateWireLengthCoef(float overflow) {
  if (overflow > 1.0) {
    wire_length_coefficient_x_ = wire_length_coefficient_y_ = 0.1;
  } else if (overflow < 0.1) {
    wire_length_coefficient_x_ = wire_length_coefficient_y_ = 10.0;
  } else {
    wire_length_coefficient_x_ = wire_length_coefficient_y_
        = 1.0 / pow(10.0, (overflow - 0.1) * 20 / 9.0 - 1.0);
  }
  wire_length_coefficient_x_ *= base_wire_length_coefficient_;
  wire_length_coefficient_y_ *= base_wire_length_coefficient_;
}
float Chip::NesterovPlacer::getPhiCoef(float scaledDiffHpwl) {
  float retCoef = (scaledDiffHpwl < 0)
                  ? max_phi_coef_
                  : max_phi_coef_
                      * pow(max_phi_coef_, scaledDiffHpwl * -1.0);
  retCoef = std::max(minPhiCoef, retCoef);
  return retCoef;
}
int64_t Chip::NesterovPlacer::getHpwl() {
  int64_t hpwl = 0;
  for (auto &gNet : net_pointers_) {
    gNet->updateBox(this->die_pointer_->getDieId());
    hpwl += gNet->hpwl();
  }
  return hpwl;
}
void Chip::NesterovPlacer::updateNextIter() {
  // swap vector pointers
  std::swap(prev_slp_coordinates_, cur_slp_coordinates_);
  std::swap(prev_slp_wire_length_grads_, cur_slp_wire_length_grads_);
  std::swap(prev_slp_density_grads_, cur_slp_density_grads_);
  std::swap(prev_slp_sum_grads_, cur_slp_sum_grads_);

  // Prevent locked instances from moving
  const auto &gCells = instance_pointers_;
  for (size_t k = 0; k < gCells.size(); ++k) {
    if (gCells[k]->isInstance() && gCells[k]->isLocked()) {
      next_slp_coordinates_[k] = cur_slp_coordinates_[k];
      next_slp_wire_length_grads_[k] = cur_slp_wire_length_grads_[k];
      next_slp_density_grads_[k] = cur_slp_density_grads_[k];
      next_slp_sum_grads_[k] = cur_slp_sum_grads_[k];

      next_coordinates_[k] = cur_coordinates_[k];
    }
  }

  std::swap(cur_slp_coordinates_, next_slp_coordinates_);
  std::swap(cur_slp_wire_length_grads_, next_slp_wire_length_grads_);
  std::swap(cur_slp_density_grads_, next_slp_density_grads_);
  std::swap(cur_slp_sum_grads_, next_slp_sum_grads_);

  std::swap(cur_coordinates_, next_coordinates_);

  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());

  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());

  updateWireLengthCoef(sum_overflow_);
  int64_t hpwl = getHpwl();

  float phiCoef = getPhiCoef(static_cast<float>(hpwl - prev_hpwl_) / referenceHpwl);

  prev_hpwl_ = hpwl;
  density_penalty_ *= phiCoef;

/*
      // for routability densityPenalty recovery
      if (rb_->numCall() == 0) {
        density_penalty_storage_.push_back(density_penalty_);
      }
*/
}
pair<float, float> Chip::NesterovPlacer::getWireLengthPreconditioner(Instance *instance) {
  // original function: getWireLengthPreconditioner
  int binding_nums = 0;
  for (Pin *pin : instance->getPins()) {
    if (pin->getNet())
      binding_nums += 1;
  }
  return pair<float, float>{static_cast<float>(binding_nums), static_cast<float>(binding_nums)};
}
pair<float, float> Chip::NesterovPlacer::getDensityPreconditioner(Instance *gCell) {
  float areaVal
      = static_cast<float>(gCell->dx()) * static_cast<float>(gCell->dy());

  return pair<float, float>{areaVal, areaVal};

}
std::pair<int, int> Chip::NesterovPlacer::getDensityMinMaxIdxX(Instance *gcell) {
  int lowerIdx = (gcell->dLx() - lx()) / bin_size_x_;
  int upperIdx = (fastModulo((gcell->dUx() - lx()), bin_size_x_) == 0)
                 ? (gcell->dUx() - lx()) / bin_size_x_
                 : (gcell->dUx() - lx()) / bin_size_x_ + 1;

  upperIdx = std::min(upperIdx, bin_cnt_x_);
  return std::make_pair(lowerIdx, upperIdx);
}
std::pair<int, int> Chip::NesterovPlacer::getDensityMinMaxIdxY(Instance *gcell) {
  int lowerIdx = (gcell->dLy() - ly()) / bin_size_y_;
  int upperIdx = (fastModulo((gcell->dUy() - ly()), bin_size_y_) == 0)
                 ? (gcell->dUy() - ly()) / bin_size_y_
                 : (gcell->dUy() - ly()) / bin_size_y_ + 1;

  upperIdx = std::min(upperIdx, bin_cnt_y_);
  return std::make_pair(lowerIdx, upperIdx);
}
float Chip::NesterovPlacer::getOverlapDensityArea(Chip::NesterovPlacer::Bin *bin, Instance *cell) {
  int rectLx = max(bin->lx(), static_cast<int>(cell->dLx()));
  int rectLy = max(bin->ly(), static_cast<int>(cell->dLy()));
  int rectUx = min(bin->ux(), static_cast<int>(cell->dUx()));
  int rectUy = min(bin->uy(), static_cast<int>(cell->dUy()));
  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  } else {
    return static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
  }
}
pair<float, float> Chip::NesterovPlacer::getDensityGradient(Instance *gCell) {
  std::pair<int, int> pairX = getDensityMinMaxIdxX(gCell);
  std::pair<int, int> pairY = getDensityMinMaxIdxY(gCell);

  pair<float, float> electroForce;

  for (int i = pairX.first; i < pairX.second; i++) {
    for (int j = pairY.first; j < pairY.second; j++) {
      Bin *bin = bins_.at(j * bin_cnt_x_ + i);
      float overlapArea
          = getOverlapDensityArea(bin, gCell) * gCell->densityScale();

      electroForce.first += overlapArea * bin->electroForceX();
      electroForce.second += overlapArea * bin->electroForceY();
    }
  }
  return electroForce;

}
void Chip::NesterovPlacer::updateGradients(vector<pair<float, float>> &sumGrads,
                                           vector<pair<float, float>> &wireLengthGrads,
                                           vector<pair<float, float>> &densityGrads) {
  wire_length_grad_sum_ = 0;
  density_grad_sum_ = 0;

  float gradSum = 0;

  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *cell = instance_pointers_.at(i);
    wireLengthGrads[i] = getWireLengthGradientWA(cell, wire_length_coefficient_x_, wire_length_coefficient_y_);
    densityGrads[i] = getDensityGradient(cell);

    // Different compiler has different results on the following formula.
    // e.g. wire_length_grad_sum_ += fabs(~~.x) + fabs(~~.y);
    //
    // To prevent instability problem,
    // I partitioned the fabs(~~.x) + fabs(~~.y) as two terms.
    //
    wire_length_grad_sum_ += fabs(wireLengthGrads[i].first);
    wire_length_grad_sum_ += fabs(wireLengthGrads[i].second);

    density_grad_sum_ += fabs(densityGrads[i].first);
    density_grad_sum_ += fabs(densityGrads[i].second);

    sumGrads[i].first = wireLengthGrads[i].first + density_penalty_ * densityGrads[i].first;
    sumGrads[i].second = wireLengthGrads[i].second + density_penalty_ * densityGrads[i].second;

    pair<float, float> wire_length_preconditioner = getWireLengthPreconditioner(cell);
    pair<float, float> density_preconditioner = getDensityPreconditioner(cell);

    pair<float, float> sum_preconditioner(
        wire_length_preconditioner.first + density_penalty_ * density_preconditioner.first,
        wire_length_preconditioner.second + density_penalty_ * density_preconditioner.second);

    if (sum_preconditioner.first <= minPreconditioner) {
      sum_preconditioner.first = minPreconditioner;
    }

    if (sum_preconditioner.second <= minPreconditioner) {
      sum_preconditioner.second = minPreconditioner;
    }

    sumGrads[i].first /= sum_preconditioner.first;
    sumGrads[i].second /= sum_preconditioner.second;

    gradSum += fabs(sumGrads[i].first) + fabs(sumGrads[i].second);
  }

  // sometimes wirelength gradient is zero when design is too small
  if (wire_length_grad_sum_ == 0
      && recursion_cnt_wl_coef_ < maxRecursionInitSLPCoef) {
    wire_length_coefficient_x_ *= 0.5;
    wire_length_coefficient_y_ *= 0.5;
    base_wire_length_coefficient_ *= 0.5;

    // update WL forces
    updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

    // recursive call again with smaller wirelength coef
    recursion_cnt_wl_coef_++;
  }

  // divergence detection on
  // Wirelength / density gradient calculation
  if (isnan(wire_length_grad_sum_) || isinf(wire_length_grad_sum_)
      || isnan(density_grad_sum_) || isinf(density_grad_sum_)) {
    is_diverged_ = true;
    diverge_msg_ = "RePlAce diverged at wire/density gradient Sum.";
    diverge_code_ = 306;
  }
}
float Chip::NesterovPlacer::getStepLength() {
  float coordiDistance = getDistance(prev_slp_coordinates_, cur_slp_coordinates_);
  float gradDistance = getDistance(prev_slp_sum_grads_, cur_slp_sum_grads_);
  return coordiDistance / gradDistance;
}
float Chip::NesterovPlacer::getDistance(vector<pair<float, float>> a, vector<pair<float, float>> b) {
  float sumDistance = 0.0f;
  for (size_t i = 0; i < a.size(); i++) {
    sumDistance += (a[i].first - b[i].first) * (a[i].first - b[i].first);
    sumDistance += (a[i].second - b[i].second) * (a[i].second - b[i].second);
  }

  return sqrt(sumDistance / (2.0 * a.size()));
}
pair<float, float> Chip::NesterovPlacer::getWireLengthGradientWA(Instance *gCell,
                                                                 float wlCoeffX,
                                                                 float wlCoeffY) const {
  pair<float, float> gradientPair;

  for (auto &gPin : gCell->getPins()) {
    if (gPin->getNet() == nullptr)
      // pass the floating pins
      continue;

    auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);

    // apply timing/custom net weight
    tmpPair.first *= gPin->getNet()->totalWeight();
    tmpPair.second *= gPin->getNet()->totalWeight();

    gradientPair.first += tmpPair.first;
    gradientPair.second += tmpPair.second;
  }

  // return sum
  assert(isnan(gradientPair.first) == false);
  assert(isnan(gradientPair.second) == false);
  return gradientPair;
}
pair<float, float> Chip::NesterovPlacer::getWireLengthGradientPinWA(Pin *gPin, float wlCoeffX, float wlCoeffY) {
  double gradientMinX = 0, gradientMinY = 0;
  double gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if (gPin->hasMinExpSumX()) {
    // from Net.
    double waExpMinSumX = gPin->getNet()->waExpMinSumX();
    double waXExpMinSumX = gPin->getNet()->waXExpMinSumX();

    gradientMinX
        = static_cast<double>(waExpMinSumX * (gPin->minExpSumX() * (1.0 - wlCoeffX * static_cast<double>(gPin->cx())))
        + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX) / (waExpMinSumX * waExpMinSumX);
  }

  // max x
  if (gPin->hasMaxExpSumX()) {
    double waExpMaxSumX = gPin->getNet()->waExpMaxSumX();
    double waXExpMaxSumX = gPin->getNet()->waXExpMaxSumX();

    gradientMaxX
        = static_cast<double>(waExpMaxSumX * (gPin->maxExpSumX() * (1.0 + wlCoeffX * static_cast<double>(gPin->cx())))
        - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX) / (waExpMaxSumX * waExpMaxSumX);
  }

  // min y
  if (gPin->hasMinExpSumY()) {
    double waExpMinSumY = gPin->getNet()->waExpMinSumY();
    double waYExpMinSumY = gPin->getNet()->waYExpMinSumY();

    gradientMinY
        = static_cast<double>(waExpMinSumY * (gPin->minExpSumY() * (1.0 - wlCoeffY * static_cast<double>(gPin->cy())))
        + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY) / (waExpMinSumY * waExpMinSumY);
  }

  // max y
  if (gPin->hasMaxExpSumY()) {
    double waExpMaxSumY = gPin->getNet()->waExpMaxSumY();
    double waYExpMaxSumY = gPin->getNet()->waYExpMaxSumY();

    gradientMaxY
        = static_cast<double>(waExpMaxSumY * (gPin->maxExpSumY() * (1.0 + wlCoeffY * static_cast<double>(gPin->cy())))
        - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY) / (waExpMaxSumY * waExpMaxSumY);
  }

  assert(!isnan(gradientMaxX));
  assert(!isnan(gradientMaxY));
  assert(!isnan(gradientMinX));
  assert(!isnan(gradientMinY));
  return pair<float, float>{gradientMinX - gradientMaxX, gradientMinY - gradientMaxY};
}
void Chip::NesterovPlacer::initSLPStepsVars() {
  const int instance_num = instance_pointers_.size();
  cur_slp_coordinates_.resize(instance_num);
  cur_slp_wire_length_grads_.resize(instance_num);
  cur_slp_density_grads_.resize(instance_num);
  cur_slp_sum_grads_.resize(instance_num);
  next_slp_coordinates_.resize(instance_num);
  next_slp_wire_length_grads_.resize(instance_num);
  next_slp_density_grads_.resize(instance_num);
  next_slp_sum_grads_.resize(instance_num);
  prev_slp_coordinates_.resize(instance_num);
  prev_slp_wire_length_grads_.resize(instance_num);
  prev_slp_density_grads_.resize(instance_num);
  prev_slp_sum_grads_.resize(instance_num);
  cur_coordinates_.resize(instance_num);
  next_coordinates_.resize(instance_num);
  init_coordinates_.resize(instance_num);
}
void Chip::NesterovPlacer::updateBinsCellDensityArea(vector<Instance *> cells) {
  // clear the Bin-area info
  for (auto &bin : bins_) {
    bin->setInstPlacedArea(0);
    bin->setInstPlacedAreaUnscaled(0);
    bin->setFillerArea(0);
  }

  for (auto &cell : cells) {
    std::pair<int, int> pairX = getDensityMinMaxIdxX(cell);
    std::pair<int, int> pairY = getDensityMinMaxIdxY(cell);

    // The following function is critical runtime hotspot
    // for global placer.
    //
    if (cell->isInstance()) {
      // macro should have
      // scale-down with target-density
      if (cell->isMacroInstance()) {
        for (int i = pairX.first; i < pairX.second; i++) {
          for (int j = pairY.first; j < pairY.second; j++) {
            Bin *bin = bins_[j * bin_cnt_x_ + i];

            const float scaledAvea = getOverlapDensityArea(bin, cell)
                * cell->densityScale()
                * bin->targetDensity();
            bin->addInstPlacedArea(scaledAvea);
            bin->addInstPlacedAreaUnscaled(scaledAvea);
          }
        }
      }
        // normal cells
      else if (cell->isStdInstance()) {
        for (int i = pairX.first; i < pairX.second; i++) {
          for (int j = pairY.first; j < pairY.second; j++) {
            Bin *bin = bins_[j * bin_cnt_x_ + i];
            const float scaledArea
                = getOverlapDensityArea(bin, cell) * cell->densityScale();
            bin->addInstPlacedArea(scaledArea);
            bin->addInstPlacedAreaUnscaled(scaledArea);
          }
        }
      }
    } else if (cell->isFiller()) {
      for (int i = pairX.first; i < pairX.second; i++) {
        for (int j = pairY.first; j < pairY.second; j++) {
          Bin *bin = bins_[j * bin_cnt_x_ + i];
          bin->addFillerArea(getOverlapDensityArea(bin, cell)
                                 * cell->densityScale());
        }
      }
    }
  }

  overflow_area_ = 0;
  overflow_area_unscaled_ = 0;
  // update density and overflowArea
  // for nesterov use and FFT library
  for (auto &bin : bins_) {
    int64_t binArea = bin->binArea();
    const float scaledBinArea = static_cast<float>(binArea * bin->targetDensity());
    bin->setDensity((static_cast<float>(bin->instPlacedArea())
        + static_cast<float>(bin->fillerArea()) + static_cast<float>(bin->nonPlaceArea())) / scaledBinArea);

    overflow_area_ += std::max(0.0f,
                               static_cast<float>(bin->instPlacedArea()) + static_cast<float>(bin->nonPlaceArea())
                                   - scaledBinArea);

    overflow_area_unscaled_ += std::max(
        0.0f,
        static_cast<float>(bin->instPlacedAreaUnscaled())
            + static_cast<float>(bin->nonPlaceAreaUnscaled()) - scaledBinArea);
  }
}
void Chip::NesterovPlacer::setDensityValuesAsDefault() {
  for (Instance *instance : instance_pointers_) {
    instance->setDensityValueAsDefault();
  }
}
void Chip::NesterovPlacer::updateDensityCoordiLayoutInside(Instance *gCell) {
  float targetLx = gCell->dLx();
  float targetLy = gCell->dLy();

  if (targetLx < lx()) {
    targetLx = lx();
  }

  if (targetLy < ly()) {
    targetLy = ly();
  }

  if (targetLx + gCell->getDensityDeltaX() > ux()) {
    targetLx = ux() - gCell->getDensityDeltaX();
  }

  if (targetLy + gCell->getDensityDeltaY() > uy()) {
    targetLy = uy() - gCell->getDensityDeltaY();
  }
  gCell->setDensityLocation(targetLx, targetLy);

}
void Chip::NesterovPlacer::updateInitialPrevSLPCoordi() {
  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *curGCell = instance_pointers_.at(i);

    float prevCoordiX
        = cur_slp_coordinates_[i].first
            - initialPrevCoordiUpdateCoef * cur_slp_sum_grads_[i].first;

    float prevCoordiY
        = cur_slp_coordinates_[i].second
            - initialPrevCoordiUpdateCoef * cur_slp_sum_grads_[i].second;

    pair<float, float> newCoordi(getDensityCoordiLayoutInsideX(curGCell, prevCoordiX),
                                 getDensityCoordiLayoutInsideY(curGCell, prevCoordiY));

    prev_slp_coordinates_[i] = newCoordi;
  }

}
void Chip::NesterovPlacer::updateDB() {
  for (Instance *instance : instance_pointers_) {
    instance->setCoordinate(instance->getDensityCenterX(), instance->getDensityCenterY());
  }
}
int Chip::NesterovPlacer::getMaxNesterovIter() const {
  return max_nesterov_iter_;
}
double fastExp(float a) {
  a = 1.0f + a / 1024.0f;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  return a;
}
} // VLSI_backend