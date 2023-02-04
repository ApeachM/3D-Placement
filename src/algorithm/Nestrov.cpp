#include <utility>

#include "Chip.h"

using namespace std;

namespace VLSI_backend {
Chip::NestrovPlacer::NestrovPlacer(odb::dbDatabase *db_database,
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

}
bool Chip::NestrovPlacer::initNestrovPlace() {
  if (!is_base_initialized) {
    // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L1054
    // bool Replace::initNesterovPlace()
    is_base_initialized = true;

    if (instance_pointers_.empty()) {
      cout << "No placeable instance - skipping placement." << endl;
      return false;
    }
    setInstancesArea();

    // void NesterovBase::init()
    // Set a fixed seed
    srand(42);
    int dbu_per_micron = db_database_->getChip()->getBlock()->getDbUnitsPerMicron();

    for (Instance *instance : instance_pointers_) {
      // For any cell, add a random noise between -1 and 1 microns to each of its
      // x and y components. This is added to make it very unlikely that identical
      // cells connected in parallel do not start at the exact same position and
      // consequently shadow each other throughout the entire placement process
      int x_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
      int y_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
      instance->setCoordinate(instance->getCoordinate().first + x_offset, instance->getCoordinate().second + y_offset);
    }

    initFillerCells();
    initBins();

    // initialize fft structure based on bins
    fft_ = new gpl::FFT(binCntX_, binCntY_, binSizeX_, binSizeY_);

    for (Instance *instance : instance_pointers_) {
      instance->setDensityValueAsDefault();
    }
    // setDensityValuesAsDefault();
    updateDensitySize();
  }

  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L138
  // void NesterovPlace::init()

  initSLPStepsVars();

  for (int idx = 0; idx < instance_pointers_.size(); ++idx) {
    Instance *gCell = instance_pointers_.at(idx);
    updateDensityCoordiLayoutInside(gCell);
    curSLPCoordi_[idx] = prevSLPCoordi_[idx] = curCoordi_[idx] = initCoordi_[idx] =
        pair<float, float>{gCell->dCx(), gCell->dCy()};
  }
  // bin
  updateGCellDensityCenterLocation(curSLPCoordi_);
  prevHpwl_ = getHpwl();
  cout << "[replace] np init: InitialHPWL: " << prevHpwl_ << endl;
  // FFT update
  updateDensityForceBin();
  baseWireLengthCoef_ = initWireLengthCoef / (static_cast<float>(binSizeX_ + binSizeY_) * 0.5);
  cout << "[replace] np init: BaseWireLengthCoef: " << baseWireLengthCoef_ << endl;
  sumOverflow_ = static_cast<float>(overflowArea_) / static_cast<float>(nesterovInstsArea());
  sumOverflowUnscaled_ = static_cast<float>(overflowAreaUnscaled_) / static_cast<float>(nesterovInstsArea());
  cout << "[replace] np init: OverflowArea: " << overflowArea_ << endl;
  cout << "[replace] np init: NesterovInstArea: " << nesterovInstsArea() << endl;
  cout << "[replace] np init: InitSumOverflow: " << sumOverflowUnscaled_ << endl;
  updateWireLengthCoef(sumOverflow_);
  cout << "[replace] np init: WireLengthCoef: " << wireLengthCoefX_ << endl;

  // WL update
  updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

  // fill in curSLPSumGrads_, curSLPWireLengthGrads_, curSLPDensityGrads_
  updateGradients(curSLPSumGrads_, curSLPWireLengthGrads_, curSLPDensityGrads_);

  if (isDiverged_) { return false; }

  // approximately fill in prevSLPCoordi_ to calculate lc vars
  updateInitialPrevSLPCoordi();

  // bin, FFT, wlen update with prevSLPCoordi.
  updateGCellDensityCenterLocation(prevSLPCoordi_);
  updateDensityForceBin();
  updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

  // update previSumGrads_, prevSLPWireLengthGrads_, prevSLPDensityGrads_
  updateGradients(prevSLPSumGrads_, prevSLPWireLengthGrads_, prevSLPDensityGrads_);

  if (isDiverged_) { return false; }

  cout << "[replace] np init: WireLengthGradSum: " << wireLengthGradSum_ << endl;
  cout << "[replace] np init: DensityGradSum: " << densityGradSum_ << endl;
  densityPenalty_ = (wireLengthGradSum_ / densityGradSum_) * initDensityPenalty;
  cout << "[replace] np init: InitDensityPenalty: " << densityPenalty_ << endl;
  sumOverflow_ = static_cast<float>(overflowArea_) / static_cast<float>(nesterovInstsArea());
  sumOverflowUnscaled_ = static_cast<float>(overflowAreaUnscaled_) / static_cast<float>(nesterovInstsArea());
  cout << "[replace] np init: PrevSumOverflow: " << sumOverflowUnscaled_ << endl;
  stepLength_ = getStepLength();
  cout << "[replace] np init: InitialStepLength: " << stepLength_ << endl;

  if ((isnan(stepLength_) || isinf(stepLength_)) && recursionCntInitSLPCoef_ < maxRecursionInitSLPCoef) {
    initialPrevCoordiUpdateCoef *= 10;
    cout << "np init: steplength = 0 detected. Rerunning Nesterov::init() with initPrevSLPCoef: "
         << initialPrevCoordiUpdateCoef << endl;
    recursionCntInitSLPCoef_++;
    initNestrovPlace();
  }

  if (isnan(stepLength_) || isinf(stepLength_)) {
    cout << "RePlAce diverged at initial iteration. Re-run with a smaller init_density_penalty value." << endl;
  }

}
int Chip::NestrovPlacer::doNestrovPlace(int start_iter) {
  // refer: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L482
  // int NesterovPlace::doNesterovPlace(int start_iter)

  // backTracking variable.
  float curA = 1.0;
  // divergence detection
  float minSumOverflow = 1e30;
  float hpwlWithMinSumOverflow = 1e30;

  // dynamic adjustment of max_phi_coef
  bool isMaxPhiCoefChanged = false;

  // snapshot saving detection
  bool isSnapshotSaved = false;

  // snapshot info
  vector<pair<float, float>> snapshotCoordi;
  vector<pair<float, float>> snapshotSLPCoordi;
  vector<pair<float, float>> snapshotSLPSumGrads;

  float snapshotA = 0;
  float snapshotDensityPenalty = 0;
  float snapshotStepLength = 0;
  float snapshotWlCoefX = 0, snapshotWlCoefY = 0;

  bool isDivergeTriedRevert = false;

  // density position initialization
  for (Instance *instance : instance_pointers_) {
    instance->setDLx(instance->getCoordinate().first);
    instance->setDLy(instance->getCoordinate().second);
    instance->setDUx(instance->getCoordinate().first + instance->getWidth());
    instance->setDUy(instance->getCoordinate().second + instance->getHeight());
  }

  // Core Nesterov Loop
  int iter = start_iter;
  for (; iter < maxNesterovIter; ++iter) {
    cout << "[replace] np: Iter: " << iter << endl;

    float prevA = curA;
    // here, prevA is a_(k), curA is a_(k+1)
    // See, the ePlace-MS paper's Algorithm 1
    curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5;
    // coeff is (a_k - 1) / ( a_(k+1) ) in paper.
    float coeff = (prevA - 1.0) / curA;

    // Back-Tracking loop
    int numBackTrak;
    for (numBackTrak = 0; numBackTrak < maxBackTrack; ++numBackTrak) {
      // fill in nextCoordinates with given stepLength_
      // here, the instance_pointers_ includes the fillers
      for (int k = 0; k < instance_pointers_.size(); ++k) {
        pair<float, float> nextCoordi(
            curSLPCoordi_[k].first + stepLength_ * curSLPSumGrads_[k].first,
            curSLPCoordi_[k].second + stepLength_ * curSLPSumGrads_[k].second);
        pair<float, float> nextSLPCoordi(
            nextCoordi.first + coeff * (nextCoordi.first - curCoordi_[k].first),
            nextCoordi.second + coeff * (nextCoordi.second - curCoordi_[k].second));
        Instance *current_instance = instance_pointers_.at(k);

        nextCoordi_.at(k) = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, nextCoordi.first),
            getDensityCoordiLayoutInsideY(current_instance, nextCoordi.second));
        nextSLPCoordi_[k] = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, nextSLPCoordi.first),
            getDensityCoordiLayoutInsideY(current_instance, nextSLPCoordi.second));
      }
      updateGCellDensityCenterLocation(nextSLPCoordi_);
      updateDensityForceBin();
      updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);
      updateGradients(nextSLPSumGrads_, nextSLPWireLengthGrads_, nextSLPDensityGrads_);

      if (isDiverged_) {
        break;
      }

      float newStepLength = getStepLength();

      if (isnan(newStepLength) || isinf(newStepLength)) {
        isDiverged_ = true;
        divergeMsg_ = "RePlAce diverged at newStepLength.";
        divergeCode_ = 305;
        break;
      }
      if (newStepLength > stepLength_ * 0.95) {
        stepLength_ = newStepLength;
        break;
      } else if (newStepLength < 0.01) {
        stepLength_ = 0.01;
        break;
      } else {
        stepLength_ = newStepLength;
      }
    }


    // dynamic adjustment for
    // better convergence with
    // large designs
    if (!isMaxPhiCoefChanged && sumOverflowUnscaled_ < 0.35f) {
      isMaxPhiCoefChanged = true;
      maxPhiCoef *= 0.99;
    }

    if (maxBackTrack == numBackTrak) {
      cout << "Backtracking limit reached so a small step will be taken" << endl;
    }

    if (isDiverged_) {
      break;
    }

    updateNextIter();

    if (iter == 0 || (iter + 1) % 10 == 0) {
      cout << iter + 1 << "\t" << sumOverflowUnscaled_ << "\t" << prevHpwl_ << endl;
      cout << "[NesterovSolve] Iter: " << iter + 1
           << "\toverflow: " << sumOverflowUnscaled_
           << "\tHPWL: " << prevHpwl_
           << endl;
    }

    if (minSumOverflow > sumOverflowUnscaled_) {
      minSumOverflow = sumOverflowUnscaled_;
      hpwlWithMinSumOverflow = prevHpwl_;
    }
    /*
        // timing driven feature
        // do reweight on timing-critical nets.
        if (npVars_.timingDrivenMode
            && tb_->isTimingNetWeightOverflow(sumOverflow_)) {
          // update db's instance location from current density coordinates
          updateDb();

          // Call resizer's estimateRC API to fill in PEX using placed locations,
          // Call sta's API to extract worst timing paths,
          // and update GNet's weights from worst timing paths.
          //
          // See timingBase.cpp in detail
          bool shouldTdProceed = tb_->updateGNetWeights(sumOverflow_);

          // problem occured
          // escape timing driven later
          if (!shouldTdProceed) {
            npVars_.timingDrivenMode = false;
          }
        }
    */
    /*
     diverge detection on
     large max_phi_cof value + large design

     1) happen overflow < 20%
     2) Hpwl is growing
    */
    if (sumOverflowUnscaled_ < 0.3f
        && sumOverflowUnscaled_ - minSumOverflow >= 0.02f
        && hpwlWithMinSumOverflow * 1.2f < prevHpwl_) {
      divergeMsg_ = "RePlAce divergence detected. ";
      divergeMsg_ += "Re-run with a smaller max_phi_cof value.";
      divergeCode_ = 307;
      isDiverged_ = true;

      // revert back to the original rb solutions
      // one more opportunity
      /*
          if (!isDivergeTriedRevert && rb_->numCall() >= 1) {
            // get back to the working rc size
            rb_->revertGCellSizeToMinRc();
      */

      // revert back the current density penality
      curCoordi_ = snapshotCoordi;
      curSLPCoordi_ = snapshotSLPCoordi;
      curSLPSumGrads_ = snapshotSLPSumGrads;
      curA = snapshotA;
      densityPenalty_ = snapshotDensityPenalty;
      stepLength_ = snapshotStepLength;
      wireLengthCoefX_ = snapshotWlCoefX;
      wireLengthCoefY_ = snapshotWlCoefY;

      updateGCellDensityCenterLocation(curCoordi_);
      updateDensityForceBin();
      updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

      isDiverged_ = false;
      divergeCode_ = 0;
      divergeMsg_ = "";
      isDivergeTriedRevert = true;

      // turn off the RD forcely
      isRoutabilityNeed_ = false;
      /*
          } else {
            // no way to revert
            break;
          }
      */
    }
    /*
        // save snapshots for routability-driven
        if (!isSnapshotSaved && npVars_.routabilityDrivenMode
            && 0.6 >= sumOverflowUnscaled_) {
          snapshotCoordi = curCoordi_;
          snapshotSLPCoordi = curSLPCoordi_;
          snapshotSLPSumGrads = curSLPSumGrads_;
          snapshotA = curA;
          snapshotDensityPenalty = densityPenalty_;
          snapshotStepLength = stepLength_;
          snapshotWlCoefX = wireLengthCoefX_;
          snapshotWlCoefY = wireLengthCoefY_;

          isSnapshotSaved = true;
          log_->report("[NesterovSolve] Snapshot saved at iter = {}", iter);
        }
    */
    /*
        // check routability using GR
        if (routabilityDrivenMode && isRoutabilityNeed_
            && routabilityCheckOverflow >= sumOverflowUnscaled_) {
          // recover the densityPenalty values
          // if further routability-driven is needed
          std::pair<bool, bool> result = rb_->routability();
          isRoutabilityNeed_ = result.first;
          bool isRevertInitNeeded = result.second;

          // if routability is needed
          if (isRoutabilityNeed_ || isRevertInitNeeded) {
            // cutFillerCoordinates();

            // revert back the current density penality
            curCoordi_ = snapshotCoordi;
            curSLPCoordi_ = snapshotSLPCoordi;
            curSLPSumGrads_ = snapshotSLPSumGrads;
            curA = snapshotA;
            densityPenalty_ = snapshotDensityPenalty;
            stepLength_ = snapshotStepLength;
            wireLengthCoefX_ = snapshotWlCoefX;
            wireLengthCoefY_ = snapshotWlCoefY;

            nb_->updateGCellDensityCenterLocation(curCoordi_);
            nb_->updateDensityForceBin();
            nb_->updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

            // reset the divergence detect conditions
            minSumOverflow = 1e30;
            hpwlWithMinSumOverflow = 1e30;
            log_->report("[NesterovSolve] Revert back to snapshot coordi");
          }
        }
    */
    // if it reached target overflow
    if (sumOverflowUnscaled_ <= targetOverflow) {
      cout << "\"[NesterovSolve] Finished with Overflow: " << sumOverflowUnscaled_ << endl;
      break;
    }
  }
  // in all case including diverge,
  // db should be updated.

  if (isDiverged_) {
    cout << "log_->error(GPL, divergeCode_, divergeMsg_);" << endl;
  }
  return iter;

}

void Chip::NestrovPlacer::setInstancesArea() {
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
void Chip::NestrovPlacer::initFillerCells() {
  // extract average dx/dy in range (10%, 90%)
  vector<int> dx_stor;
  vector<int> dy_stor;
  dx_stor.reserve(instance_pointers_.size());
  dy_stor.reserve(instance_pointers_.size());
  for (Instance *instance : instance_pointers_) {
    dx_stor.push_back(static_cast<int>(instance->getWidth()));
    dy_stor.push_back(static_cast<int>(instance->getHeight()));
  }

  // sort
  std::sort(dx_stor.begin(), dx_stor.end());
  std::sort(dy_stor.begin(), dy_stor.end());

  // average from (10 - 90%)
  int64_t dx_sum = 0, dy_sum = 0;

  int min_idx = static_cast<int>(static_cast<float>(dx_stor.size()) * 0.05);
  int max_idx = static_cast<int>(static_cast<float>(dx_stor.size()) * 0.95);

  // when #instances are too small,
  // extracts average values in whole ranges.
  if (min_idx == max_idx) {
    min_idx = 0;
    max_idx = static_cast<int>(dx_stor.size());
  }

  for (int i = min_idx; i < max_idx; i++) {
    dx_sum += dx_stor[i];
    dy_sum += dy_stor[i];
  }

  // the avgDx and avgDy will be used as filler cells
  // width and height
  filler_width_ = static_cast<int>(dx_sum / (max_idx - min_idx));
  filler_height_ = static_cast<int>(dy_sum / (max_idx - min_idx));

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
void Chip::NestrovPlacer::initBins() {
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
  binCntX_ = foundBinCnt;
  binCntY_ = foundBinCnt;
  // log_->info(GPL, 28, "BinCnt: {} {}", binCntX_, binCntY_);
  binSizeX_ = ceil(static_cast<float>((ux_ - lx_)) / binCntX_);
  binSizeY_ = ceil(static_cast<float>((uy_ - ly_)) / binCntY_);
  // log_->info(GPL, 29, "BinSize: {} {}", binSizeX_, binSizeY_);
  binStor_.resize(binCntX_ * binCntY_);
  bins_.reserve(binCntX_ * binCntY_);
  int x = lx_, y = ly_;
  int idxX = 0, idxY = 0;
  for (auto &bin : binStor_) {
    int sizeX = (x + binSizeX_ > ux_) ? ux_ - x : binSizeX_;
    int sizeY = (y + binSizeY_ > uy_) ? uy_ - y : binSizeY_;
    bin = Bin(idxX, idxY, x, y, x + sizeX, y + sizeY, target_density_);

    // move x, y coordinates.
    x += binSizeX_;
    idxX += 1;
    if (x >= ux_) {
      y += binSizeY_;
      x = lx_;
      idxY++;
      idxX = 0;
    }
    bins_.push_back(&bin);
  }


  // only initialized once
  updateBinsNonPlaceArea();
}
void Chip::NestrovPlacer::updateBinsNonPlaceArea() {
  for (auto &bin : bins_) {
    bin->setNonPlaceArea(0);
    bin->setNonPlaceAreaUnscaled(0);
  }

  for (auto &inst : nonPlaceInsts_) {
    std::pair<int, int> pairX = getMinMaxIdxX(inst);
    std::pair<int, int> pairY = getMinMaxIdxY(inst);
    for (int i = pairX.first; i < pairX.second; i++) {
      for (int j = pairY.first; j < pairY.second; j++) {
        Bin *bin = bins_[j * binCntX_ + i];

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
float Chip::NestrovPlacer::getDensityCoordiLayoutInsideX(Instance *instance, float cx) {
  float adjVal = cx;  // adjusted value
  // will change base on each assigned binGrids.
  if (cx - instance->dDx() / 2 < die_pointer_->getLowerLeftX()) {
    adjVal = die_pointer_->getLowerLeftX() + instance->dDx() / 2;
  }
  if (cx + instance->dDx() / 2 > die_pointer_->getUpperRightX()) {
    adjVal = die_pointer_->getUpperRightX() - instance->dDx() / 2;
  }
  return adjVal;
}
float Chip::NestrovPlacer::getDensityCoordiLayoutInsideY(Instance *instance, float cy) {
  float adjVal = cy;  // adjusted value
  // will change base on each assigned binGrids.
  if (cy - instance->dDy() / 2 < die_pointer_->getLowerLeftY()) {
    adjVal = die_pointer_->getLowerLeftY() + instance->dDy() / 2;
  }
  if (cy + instance->dDy() / 2 > die_pointer_->getUpperRightY()) {
    adjVal = die_pointer_->getUpperRightY() - instance->dDy() / 2;
  }
  return adjVal;
}
void Chip::NestrovPlacer::updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates) {
  for (int idx = 0; idx < coordinates.size(); ++idx) {
    pair<float, float> coordinate = coordinates.at(idx);
    instance_pointers_.at(idx)->setDensityCenterLocation(coordinate.first, coordinate.second);
  }
  updateBinsGCellDensityArea(instance_pointers_);
}
std::pair<int, int> Chip::NestrovPlacer::getMinMaxIdxX(Instance *inst) const {
  int lowerIdx = (inst->ly() - die_pointer_->getLowerLeftY()) / binSizeY_;
  int upperIdx = (fastModulo((inst->uy() - lx()), binSizeY_) == 0)
                 ? (inst->uy() - ly()) / binSizeY_
                 : (inst->uy() - ly()) / binSizeY_ + 1;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
}
std::pair<int, int> Chip::NestrovPlacer::getMinMaxIdxY(Instance *inst) const {
  int lowerIdx = (inst->ly() - ly()) / binSizeY_;
  int upperIdx = (fastModulo((inst->uy() - ly()), binSizeY_) == 0)
                 ? (inst->uy() - ly()) / binSizeY_
                 : (inst->uy() - ly()) / binSizeY_ + 1;

  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, binCntY_));
}
int Chip::NestrovPlacer::fastModulo(const int input, const int ceil) { return input >= ceil ? input % ceil : input; }
int64_t Chip::NestrovPlacer::getOverlapArea(const Chip::NestrovPlacer::Bin *bin, Instance *inst, int dbu_per_micron) {
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
float Chip::NestrovPlacer::calculateBiVariateNormalCDF(Chip::NestrovPlacer::biNormalParameters i) {
  /*!\brief
  A function that does 2D integration to the density function of a
  bivariate normal distribution with 0 correlation.
  Essentially, the function being integrated is the product
  of 2 1D probability density functions (for x and y). The means and standard
  deviation of the probablity density functions are parametarized. In this
  function, I am using the closed-form solution of the integration. The limits
  of integration are lx->ux and ly->uy For reference: the equation that is
  being integrated is:
       (1/(2*pi*sigmaX*sigmaY))*e^(-(y-meanY)^2/(2*sigmaY*sigmaY))*e^(-(x-meanX)^2/(2*sigmaX*sigmaX))
  */
  const float x1 = (i.meanX - i.lx) / (sqrt(2) * i.sigmaX);
  const float x2 = (i.meanX - i.ux) / (sqrt(2) * i.sigmaX);

  const float y1 = (i.meanY - i.ly) / (sqrt(2) * i.sigmaY);
  const float y2 = (i.meanY - i.uy) / (sqrt(2) * i.sigmaY);

  return 0.25
      * (erf(x1) * erf(y1) + erf(x2) * erf(y2) - erf(x1) * erf(y2)
          - erf(x2) * erf(y1));
}
int64_t Chip::NestrovPlacer::getOverlapAreaUnscaled(const Chip::NestrovPlacer::Bin *bin, Instance *inst) {
  int rectLx = max(bin->lx(), inst->lx()), rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), rectUy = min(bin->uy(), inst->uy());

  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  } else {
    return static_cast<int64_t>(rectUx - rectLx)
        * static_cast<int64_t>(rectUy - rectLy);
  }
}
void Chip::NestrovPlacer::updateDensitySize() {
  for (Instance *instance : instance_pointers_) {
    float scaleX = 0, scaleY = 0;
    float densitySizeX = 0, densitySizeY = 0;
    if (instance->dx() < REPLACE_SQRT2 * binSizeX_) {
      scaleX = static_cast<float>(instance->dx())
          / static_cast<float>(REPLACE_SQRT2 * binSizeX_);
      densitySizeX = REPLACE_SQRT2 * static_cast<float>(binSizeX_);
    } else {
      scaleX = 1.0;
      densitySizeX = instance->dx();
    }

    if (instance->dy() < REPLACE_SQRT2 * binSizeY_) {
      scaleY = static_cast<float>(instance->dy())
          / static_cast<float>(REPLACE_SQRT2 * binSizeY_);
      densitySizeY = REPLACE_SQRT2 * static_cast<float>(binSizeY_);
    } else {
      scaleY = 1.0;
      densitySizeY = instance->dy();
    }

    instance->setDensitySize(densitySizeX, densitySizeY);
    instance->setDensityScale(scaleX * scaleY);
  }
}
void Chip::NestrovPlacer::updateDensityForceBin() {
  // copy density to utilize FFT
  for (auto &bin : bins_) {
    fft_->updateDensity(bin->x(), bin->y(), bin->density());
  }

  // do FFT
  fft_->doFFT();

  // update electroPhi and electroForce
  // update sumPhi_ for nesterov loop
  sumPhi_ = 0;
  for (auto &bin : bins_) {
    auto eForcePair = fft_->getElectroForce(bin->x(), bin->y());
    bin->setElectroForce(eForcePair.first, eForcePair.second);

    float electroPhi = fft_->getElectroPhi(bin->x(), bin->y());
    bin->setElectroPhi(electroPhi);

    sumPhi_ += electroPhi * static_cast<float>(bin->nonPlaceArea() + bin->instPlacedArea() + bin->fillerArea());
  }

}
void Chip::NestrovPlacer::updateWireLengthForceWA(float wlCoeffX, float wlCoeffY) {
  // clear all WA variables.
  for (Net *gNet : net_pointers_) {
    gNet->clearWaVars();
  }
  for (auto &gPin : pin_pointers_) {
    gPin->clearWaVars();
  }

  for (Net *&gNet : net_pointers_) {
    gNet->updateBox();

    for (Pin *gPin : gNet->getConnectedPins()) {
      // The WA terms are shift invariant:
      //
      //   Sum(x_i * exp(x_i))    Sum(x_i * exp(x_i - C))
      //   -----------------    = -----------------
      //   Sum(exp(x_i))          Sum(exp(x_i - C))
      //
      // So we shift to keep the exponential from overflowing
      float expMinX = static_cast<float>(gNet->lx() - gPin->cx()) * wlCoeffX;
      float expMaxX = static_cast<float>(gPin->cx() - gNet->ux()) * wlCoeffX;
      float expMinY = static_cast<float>(gNet->ly() - gPin->cy()) * wlCoeffY;
      float expMaxY = static_cast<float>(gPin->cy() - gNet->uy()) * wlCoeffY;

      // min x
      if (expMinX > minWireLengthForceBar) {
        gPin->setMinExpSumX(fastExp(expMinX));
        gNet->addWaExpMinSumX(gPin->minExpSumX());
        gNet->addWaXExpMinSumX(gPin->cx() * gPin->minExpSumX());
      }

      // max x
      if (expMaxX > minWireLengthForceBar) {
        gPin->setMaxExpSumX(fastExp(expMaxX));
        gNet->addWaExpMaxSumX(gPin->maxExpSumX());
        gNet->addWaXExpMaxSumX(gPin->cx() * gPin->maxExpSumX());
      }

      // min y
      if (expMinY > minWireLengthForceBar) {
        gPin->setMinExpSumY(fastExp(expMinY));
        gNet->addWaExpMinSumY(gPin->minExpSumY());
        gNet->addWaYExpMinSumY(gPin->cy() * gPin->minExpSumY());
      }

      // max y
      if (expMaxY > minWireLengthForceBar) {
        gPin->setMaxExpSumY(fastExp(expMaxY));
        gNet->addWaExpMaxSumY(gPin->maxExpSumY());
        gNet->addWaYExpMaxSumY(gPin->cy() * gPin->maxExpSumY());
      }
    }
  }

}
float Chip::NestrovPlacer::nesterovInstsArea() const {
  return std_instances_area_ + static_cast<int64_t>(round(macro_instances_area_ * target_density_));
}
void Chip::NestrovPlacer::updateWireLengthCoef(float overflow) {
  if (overflow > 1.0) {
    wireLengthCoefX_ = wireLengthCoefY_ = 0.1;
  } else if (overflow < 0.1) {
    wireLengthCoefX_ = wireLengthCoefY_ = 10.0;
  } else {
    wireLengthCoefX_ = wireLengthCoefY_
        = 1.0 / pow(10.0, (overflow - 0.1) * 20 / 9.0 - 1.0);
  }
  wireLengthCoefX_ *= baseWireLengthCoef_;
  wireLengthCoefY_ *= baseWireLengthCoef_;
}
float Chip::NestrovPlacer::getPhiCoef(float scaledDiffHpwl) {
  float retCoef = (scaledDiffHpwl < 0)
                  ? maxPhiCoef
                  : maxPhiCoef
                      * pow(maxPhiCoef, scaledDiffHpwl * -1.0);
  retCoef = std::max(minPhiCoef, retCoef);
  return retCoef;
}
int64_t Chip::NestrovPlacer::getHpwl() {
  int64_t hpwl = 0;
  for (auto &gNet : net_pointers_) {
    gNet->updateBox();
    hpwl += gNet->hpwl();
  }
  return hpwl;
}
void Chip::NestrovPlacer::updateNextIter() {
  // swap vector pointers
  std::swap(prevSLPCoordi_, curSLPCoordi_);
  std::swap(prevSLPWireLengthGrads_, curSLPWireLengthGrads_);
  std::swap(prevSLPDensityGrads_, curSLPDensityGrads_);
  std::swap(prevSLPSumGrads_, curSLPSumGrads_);

  // Prevent locked instances from moving
  const auto &gCells = instance_pointers_;
  for (size_t k = 0; k < gCells.size(); ++k) {
    if (gCells[k]->isInstance() && gCells[k]->isLocked()) {
      nextSLPCoordi_[k] = curSLPCoordi_[k];
      nextSLPWireLengthGrads_[k] = curSLPWireLengthGrads_[k];
      nextSLPDensityGrads_[k] = curSLPDensityGrads_[k];
      nextSLPSumGrads_[k] = curSLPSumGrads_[k];

      nextCoordi_[k] = curCoordi_[k];
    }
  }

  std::swap(curSLPCoordi_, nextSLPCoordi_);
  std::swap(curSLPWireLengthGrads_, nextSLPWireLengthGrads_);
  std::swap(curSLPDensityGrads_, nextSLPDensityGrads_);
  std::swap(curSLPSumGrads_, nextSLPSumGrads_);

  std::swap(curCoordi_, nextCoordi_);

  sumOverflow_ = static_cast<float>(overflowArea_) / static_cast<float>(nesterovInstsArea());

  sumOverflowUnscaled_ = static_cast<float>(overflowAreaUnscaled_) / static_cast<float>(nesterovInstsArea());

  updateWireLengthCoef(sumOverflow_);
  int64_t hpwl = getHpwl();

  float phiCoef = getPhiCoef(static_cast<float>(hpwl - prevHpwl_) / referenceHpwl);

  prevHpwl_ = hpwl;
  densityPenalty_ *= phiCoef;

/*
      // for routability densityPenalty recovery
      if (rb_->numCall() == 0) {
        densityPenaltyStor_.push_back(densityPenalty_);
      }
*/
}
pair<float, float> Chip::NestrovPlacer::getWireLengthPreconditioner(Instance *instance) {
  // original function: getWireLengthPreconditioner
  int binding_nums = 0;
  for (Pin* pin: instance->getPins()) {
    if (pin->getNet())
      binding_nums += 1;
  }
  return pair<float, float>{static_cast<float>(binding_nums), static_cast<float>(binding_nums)};
}
pair<float, float> Chip::NestrovPlacer::getDensityPreconditioner(Instance *gCell) {
  float areaVal
      = static_cast<float>(gCell->dx()) * static_cast<float>(gCell->dy());

  return pair<float, float>{areaVal, areaVal};

}
std::pair<int, int> Chip::NestrovPlacer::getDensityMinMaxIdxX(Instance *gcell) {
  int lowerIdx = (gcell->dLx() - lx()) / binSizeX_;
  int upperIdx = (fastModulo((gcell->dUx() - lx()), binSizeX_) == 0)
                 ? (gcell->dUx() - lx()) / binSizeX_
                 : (gcell->dUx() - lx()) / binSizeX_ + 1;

  upperIdx = std::min(upperIdx, binCntX_);
  return std::make_pair(lowerIdx, upperIdx);
}
std::pair<int, int> Chip::NestrovPlacer::getDensityMinMaxIdxY(Instance *gcell) {
  int lowerIdx = (gcell->dLy() - ly()) / binSizeY_;
  int upperIdx = (fastModulo((gcell->dUy() - ly()), binSizeY_) == 0)
                 ? (gcell->dUy() - ly()) / binSizeY_
                 : (gcell->dUy() - ly()) / binSizeY_ + 1;

  upperIdx = std::min(upperIdx, binCntY_);
  return std::make_pair(lowerIdx, upperIdx);
}
float Chip::NestrovPlacer::getOverlapDensityArea(Chip::NestrovPlacer::Bin *bin, Instance *cell) {
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
pair<float, float> Chip::NestrovPlacer::getDensityGradient(Instance *gCell) {
  std::pair<int, int> pairX = getDensityMinMaxIdxX(gCell);
  std::pair<int, int> pairY = getDensityMinMaxIdxY(gCell);

  pair<float, float> electroForce;

  for (int i = pairX.first; i < pairX.second; i++) {
    for (int j = pairY.first; j < pairY.second; j++) {
      Bin *bin = bins_.at(j * binCntX_ + i);
      float overlapArea
          = getOverlapDensityArea(bin, gCell) * gCell->densityScale();

      electroForce.first += overlapArea * bin->electroForceX();
      electroForce.second += overlapArea * bin->electroForceY();
    }
  }
  return electroForce;

}
void Chip::NestrovPlacer::updateGradients(vector<pair<float, float>> &sumGrads,
                                          vector<pair<float, float>> &wireLengthGrads,
                                          vector<pair<float, float>> &densityGrads) {
  wireLengthGradSum_ = 0;
  densityGradSum_ = 0;

  float gradSum = 0;

  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *gCell = instance_pointers_.at(i);
    wireLengthGrads[i] = getWireLengthGradientWA(gCell, wireLengthCoefX_, wireLengthCoefY_);
    densityGrads[i] = getDensityGradient(gCell);

    // Different compiler has different results on the following formula.
    // e.g. wireLengthGradSum_ += fabs(~~.x) + fabs(~~.y);
    //
    // To prevent instability problem,
    // I partitioned the fabs(~~.x) + fabs(~~.y) as two terms.
    //
    wireLengthGradSum_ += fabs(wireLengthGrads[i].first);
    wireLengthGradSum_ += fabs(wireLengthGrads[i].second);

    densityGradSum_ += fabs(densityGrads[i].first);
    densityGradSum_ += fabs(densityGrads[i].second);

    sumGrads[i].first = wireLengthGrads[i].first + densityPenalty_ * densityGrads[i].first;
    sumGrads[i].second = wireLengthGrads[i].second + densityPenalty_ * densityGrads[i].second;

    pair<float, float> wireLengthPreCondi = getWireLengthPreconditioner(gCell);
    pair<float, float> densityPrecondi = getDensityPreconditioner(gCell);

    pair<float, float> sumPrecondi(
        wireLengthPreCondi.first + densityPenalty_ * densityPrecondi.first,
        wireLengthPreCondi.second + densityPenalty_ * densityPrecondi.second);

    if (sumPrecondi.first <= minPreconditioner) {
      sumPrecondi.first = minPreconditioner;
    }

    if (sumPrecondi.second <= minPreconditioner) {
      sumPrecondi.second = minPreconditioner;
    }

    sumGrads[i].first /= sumPrecondi.first;
    sumGrads[i].second /= sumPrecondi.second;

    gradSum += fabs(sumGrads[i].first) + fabs(sumGrads[i].second);
  }

  // sometimes wirelength gradient is zero when design is too small
  if (wireLengthGradSum_ == 0
      && recursionCntWlCoef_ < maxRecursionInitSLPCoef) {
    wireLengthCoefX_ *= 0.5;
    wireLengthCoefY_ *= 0.5;
    baseWireLengthCoef_ *= 0.5;

    // update WL forces
    updateWireLengthForceWA(wireLengthCoefX_, wireLengthCoefY_);

    // recursive call again with smaller wirelength coef
    recursionCntWlCoef_++;
  }

  // divergence detection on
  // Wirelength / density gradient calculation
  if (isnan(wireLengthGradSum_) || isinf(wireLengthGradSum_)
      || isnan(densityGradSum_) || isinf(densityGradSum_)) {
    isDiverged_ = true;
    divergeMsg_ = "RePlAce diverged at wire/density gradient Sum.";
    divergeCode_ = 306;
  }
}
float Chip::NestrovPlacer::getStepLength() {
  float coordiDistance = getDistance(prevSLPCoordi_, curSLPCoordi_);
  float gradDistance = getDistance(prevSLPSumGrads_, curSLPSumGrads_);
  return coordiDistance / gradDistance;
}
float Chip::NestrovPlacer::getDistance(vector<pair<float, float>> a, vector<pair<float, float>> b) {
  float sumDistance = 0.0f;
  for (size_t i = 0; i < a.size(); i++) {
    sumDistance += (a[i].first - b[i].first) * (a[i].first - b[i].first);
    sumDistance += (a[i].second - b[i].second) * (a[i].second - b[i].second);
  }

  return sqrt(sumDistance / (2.0 * a.size()));
}
pair<float, float> Chip::NestrovPlacer::getWireLengthGradientWA(Instance *gCell, float wlCoeffX, float wlCoeffY) const {
  pair<float, float> gradientPair;

  for (auto &gPin : gCell->getPins()) {
    if (gPin->getNet() == nullptr) {
      continue;
      // TODO
      //  should check this is right or not
    }

    auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);

    // apply timing/custom net weight
    tmpPair.first *= gPin->getNet()->totalWeight();
    tmpPair.second *= gPin->getNet()->totalWeight();

    gradientPair.first += tmpPair.first;
    gradientPair.second += tmpPair.second;
  }

  // return sum
  return gradientPair;
}
pair<float, float> Chip::NestrovPlacer::getWireLengthGradientPinWA(Pin *gPin, float wlCoeffX, float wlCoeffY) {
  float gradientMinX = 0, gradientMinY = 0;
  float gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if (gPin->hasMinExpSumX()) {
    // from Net.
    float waExpMinSumX = gPin->getNet()->waExpMinSumX();
    float waXExpMinSumX = gPin->getNet()->waXExpMinSumX();

    gradientMinX
        = static_cast<float>(waExpMinSumX * (gPin->minExpSumX() * (1.0 - wlCoeffX * static_cast<float>(gPin->cx())))
        + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX) / (waExpMinSumX * waExpMinSumX);
  }

  // max x
  if (gPin->hasMaxExpSumX()) {
    float waExpMaxSumX = gPin->getNet()->waExpMaxSumX();
    float waXExpMaxSumX = gPin->getNet()->waXExpMaxSumX();

    gradientMaxX
        = static_cast<float>(waExpMaxSumX * (gPin->maxExpSumX() * (1.0 + wlCoeffX * static_cast<float>(gPin->cx())))
        - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX) / (waExpMaxSumX * waExpMaxSumX);
  }

  // min y
  if (gPin->hasMinExpSumY()) {
    float waExpMinSumY = gPin->getNet()->waExpMinSumY();
    float waYExpMinSumY = gPin->getNet()->waYExpMinSumY();

    gradientMinY
        = static_cast<float>(waExpMinSumY * (gPin->minExpSumY() * (1.0 - wlCoeffY * static_cast<float>(gPin->cy())))
        + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY) / (waExpMinSumY * waExpMinSumY);
  }

  // max y
  if (gPin->hasMaxExpSumY()) {
    float waExpMaxSumY = gPin->getNet()->waExpMaxSumY();
    float waYExpMaxSumY = gPin->getNet()->waYExpMaxSumY();

    gradientMaxY
        = static_cast<float>(waExpMaxSumY * (gPin->maxExpSumY() * (1.0 + wlCoeffY * static_cast<float>(gPin->cy())))
        - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY) / (waExpMaxSumY * waExpMaxSumY);
  }

  return pair<float, float>{gradientMinX - gradientMaxX, gradientMinY - gradientMaxY};
}
void Chip::NestrovPlacer::initSLPStepsVars() {
  const int instance_num = instance_pointers_.size();
  curSLPCoordi_.resize(instance_num);
  curSLPWireLengthGrads_.resize(instance_num);
  curSLPDensityGrads_.resize(instance_num);
  curSLPSumGrads_.resize(instance_num);
  nextSLPCoordi_.resize(instance_num);
  nextSLPWireLengthGrads_.resize(instance_num);
  nextSLPDensityGrads_.resize(instance_num);
  nextSLPSumGrads_.resize(instance_num);
  prevSLPCoordi_.resize(instance_num);
  prevSLPWireLengthGrads_.resize(instance_num);
  prevSLPDensityGrads_.resize(instance_num);
  prevSLPSumGrads_.resize(instance_num);
  curCoordi_.resize(instance_num);
  nextCoordi_.resize(instance_num);
  initCoordi_.resize(instance_num);
}
void Chip::NestrovPlacer::updateBinsGCellDensityArea(vector<Instance *> cells) {
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
            Bin *bin = bins_[j * binCntX_ + i];

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
            Bin *bin = bins_[j * binCntX_ + i];
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
          Bin *bin = bins_[j * binCntX_ + i];
          bin->addFillerArea(getOverlapDensityArea(bin, cell)
                                 * cell->densityScale());
        }
      }
    }
  }

  overflowArea_ = 0;
  overflowAreaUnscaled_ = 0;
  // update density and overflowArea
  // for nesterov use and FFT library
  for (auto &bin : bins_) {
    int64_t binArea = bin->binArea();
    const float scaledBinArea
        = static_cast<float>(binArea * bin->targetDensity());
    bin->setDensity((static_cast<float>(bin->instPlacedArea())
        + static_cast<float>(bin->fillerArea())
        + static_cast<float>(bin->nonPlaceArea()))
                        / scaledBinArea);

    overflowArea_ += std::max(0.0f,
                              static_cast<float>(bin->instPlacedArea())
                                  + static_cast<float>(bin->nonPlaceArea())
                                  - scaledBinArea);

    overflowAreaUnscaled_ += std::max(
        0.0f,
        static_cast<float>(bin->instPlacedAreaUnscaled())
            + static_cast<float>(bin->nonPlaceAreaUnscaled()) - scaledBinArea);
  }
}
void Chip::NestrovPlacer::setDensityValuesAsDefault() {
  for (Instance *instance : instance_pointers_) {
    instance->setDensityValueAsDefault();
  }
}
void Chip::NestrovPlacer::updateDensityCoordiLayoutInside(Instance *gCell) {
  float targetLx = gCell->dLx();
  float targetLy = gCell->dLy();

  if (targetLx < lx()) {
    targetLx = lx();
  }

  if (targetLy < ly()) {
    targetLy = ly();
  }

  if (targetLx + gCell->dDx() > ux()) {
    targetLx = ux() - gCell->dDx();
  }

  if (targetLy + gCell->dDy() > uy()) {
    targetLy = uy() - gCell->dDy();
  }
  gCell->setDensityLocation(targetLx, targetLy);

}
void Chip::NestrovPlacer::updateInitialPrevSLPCoordi() {
  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *curGCell = instance_pointers_.at(i);

    float prevCoordiX
        = curSLPCoordi_[i].first
            - initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].first;

    float prevCoordiY
        = curSLPCoordi_[i].second
            - initialPrevCoordiUpdateCoef * curSLPSumGrads_[i].second;

    pair<float, float> newCoordi(getDensityCoordiLayoutInsideX(curGCell, prevCoordiX),
                                 getDensityCoordiLayoutInsideY(curGCell, prevCoordiY));

    prevSLPCoordi_[i] = newCoordi;
  }

}
float fastExp(float a) {
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