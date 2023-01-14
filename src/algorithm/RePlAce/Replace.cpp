#include "Replace.h"
#include "nesterovBase.h"

namespace gpl {

Replace::Replace()
    : db_(nullptr),
      pb_(nullptr),
      nb_(nullptr),
      ip_(nullptr),
      np_(nullptr),
      initialPlaceMaxIter_(20),
      initialPlaceMinDiffLength_(1500),
      initialPlaceMaxSolverIter_(100),
      initialPlaceMaxFanout_(200),
      initialPlaceNetWeightScale_(800),
      forceCPU_(false),
      nesterovPlaceMaxIter_(5000),
      binGridCntX_(0),
      binGridCntY_(0),
      overflow_(0.1),
      density_(1.0),
      initDensityPenalityFactor_(0.00008),
      initWireLengthCoef_(0.25),
      minPhiCoef_(0.95),
      maxPhiCoef_(1.05),
      referenceHpwl_(446000000),
      routabilityCheckOverflow_(0.20),
      routabilityMaxDensity_(0.99),
      routabilityTargetRcMetric_(1.25),
      routabilityInflationRatioCoef_(2.5),
      routabilityMaxInflationRatio_(2.5),
      routabilityRcK1_(1.0),
      routabilityRcK2_(1.0),
      routabilityRcK3_(0.0),
      routabilityRcK4_(0.0),
      routabilityMaxBloatIter_(1),
      routabilityMaxInflationIter_(4),
      timingNetWeightMax_(1.9),
      timingDrivenMode_(true),
      routabilityDrivenMode_(true),
      uniformTargetDensityMode_(false),
      skipIoMode_(false),
      padLeft_(0),
      padRight_(0) {};

Replace::~Replace() {

}
void Replace::init() {

}
void Replace::reset() {

  // two pointers should not be freed.
  db_ = nullptr;

  ip_.reset();
  np_.reset();

  pb_.reset();
  nb_.reset();

  initialPlaceMaxIter_ = 20;
  initialPlaceMinDiffLength_ = 1500;
  initialPlaceMaxSolverIter_ = 100;
  initialPlaceMaxFanout_ = 200;
  initialPlaceNetWeightScale_ = 800;
  forceCPU_ = false;

  nesterovPlaceMaxIter_ = 5000;
  binGridCntX_ = binGridCntY_ = 0;
  overflow_ = 0.1;
  density_ = 1.0;
  initDensityPenalityFactor_ = 0.00008;
  initWireLengthCoef_ = 0.25;
  minPhiCoef_ = 0.95;
  maxPhiCoef_ = 1.05;
  referenceHpwl_ = 446000000;

  routabilityCheckOverflow_ = 0.20;
  routabilityMaxDensity_ = 0.99;
  routabilityTargetRcMetric_ = 1.25;
  routabilityInflationRatioCoef_ = 2.5;
  routabilityMaxInflationRatio_ = 2.5;
  routabilityRcK1_ = routabilityRcK2_ = 1.0;
  routabilityRcK3_ = routabilityRcK4_ = 0.0;
  routabilityMaxBloatIter_ = 1;
  routabilityMaxInflationIter_ = 4;

  timingDrivenMode_ = true;
  routabilityDrivenMode_ = true;
  uniformTargetDensityMode_ = false;
  skipIoMode_ = false;

  padLeft_ = padRight_ = 0;

  timingNetWeightOverflows_.clear();
  timingNetWeightOverflows_.shrink_to_fit();
  timingNetWeightMax_ = 1.9;

}
void Replace::setDb(odb::dbDatabase *odb) {
  db_ = odb;
}
void Replace::doIncrementalPlace() {
  PlacerBaseVars pbVars;
  pbVars.padLeft = padLeft_;
  pbVars.padRight = padRight_;
  pbVars.skipIoMode = skipIoMode_;

  pb_ = std::make_shared<PlacerBase>(db_, pbVars, log_);

  // Lock down already placed objects
  int locked_cnt = 0;
  int unplaced_cnt = 0;
  auto block = db_->getChip()->getBlock();
  for (auto inst : block->getInsts()) {
    auto status = inst->getPlacementStatus();
    if (status == odb::dbPlacementStatus::PLACED) {
      pb_->dbToPb(inst)->lock();
      ++locked_cnt;
    } else if (!status.isPlaced()) {
      ++unplaced_cnt;
    }
  }

  if (unplaced_cnt == 0) {
    // Everything was already placed so we do the old incremental mode
    // which just skips initial placement and runs nesterov.
    pb_->unlockAll();
    doNesterovPlace();
    return;
  }


  // Roughly place the unplaced objects (allow more overflow).
  // Limit iterations to prevent objects drifting too far or
  // non-convergence.
  constexpr float rough_oveflow = 0.2f;
  float previous_overflow = overflow_;
  setTargetOverflow(std::max(rough_oveflow, overflow_));
  doInitialPlace();

  int previous_max_iter = nesterovPlaceMaxIter_;
  initNesterovPlace();
  setNesterovPlaceMaxIter(300);
  int iter = doNesterovPlace();
  setNesterovPlaceMaxIter(previous_max_iter);

  // Finish the overflow resolution from the rough placement
  pb_->unlockAll();

  setTargetOverflow(previous_overflow);
  if (previous_overflow < rough_oveflow) {
    doNesterovPlace(iter + 1);
  }

}
void Replace::doInitialPlace() {
  if (pb_ == nullptr) {
    PlacerBaseVars pbVars;
    pbVars.padLeft = padLeft_;
    pbVars.padRight = padRight_;
    pbVars.skipIoMode = skipIoMode_;

    pb_ = std::make_shared<PlacerBase>(db_, pbVars, log_);
  }

  InitialPlaceVars ipVars;
  ipVars.maxIter = initialPlaceMaxIter_;
  ipVars.minDiffLength = initialPlaceMinDiffLength_;
  ipVars.maxSolverIter = initialPlaceMaxSolverIter_;
  ipVars.maxFanout = initialPlaceMaxFanout_;
  ipVars.netWeightScale = initialPlaceNetWeightScale_;
  ipVars.forceCPU = forceCPU_;

  std::unique_ptr<InitialPlace> ip(new InitialPlace(ipVars, pb_, log_));
  ip_ = std::move(ip);
  ip_->doBicgstabPlace();
}
bool Replace::initNesterovPlace() {
  if (!pb_) {
    PlacerBaseVars pbVars;
    pbVars.padLeft = padLeft_;
    pbVars.padRight = padRight_;
    pbVars.skipIoMode = skipIoMode_;

    pb_ = std::make_shared<PlacerBase>(db_, pbVars, log_);
  }

  if (pb_->placeInsts().size() == 0) {
/*
    log_->warn(GPL, 136, "No placeable instances - skipping placement.");
*/
    std::cout << "No placeable instances - skipping placement." << std::endl;
    return false;
  }

  if (!nb_) {
    NesterovBaseVars nbVars;
    nbVars.targetDensity = density_;

    if (binGridCntX_ != 0) {
      nbVars.isSetBinCntX = 1;
      nbVars.binCntX = binGridCntX_;
    }

    if (binGridCntY_ != 0) {
      nbVars.isSetBinCntY = 1;
      nbVars.binCntY = binGridCntY_;
    }

    nbVars.useUniformTargetDensity = uniformTargetDensityMode_;

    nb_ = std::make_shared<NesterovBase>(nbVars, pb_, log_);
  }

  if (!np_) {
    NesterovPlaceVars npVars;

    npVars.minPhiCoef = minPhiCoef_;
    npVars.maxPhiCoef = maxPhiCoef_;
    npVars.referenceHpwl = referenceHpwl_;
    npVars.routabilityCheckOverflow = routabilityCheckOverflow_;
    npVars.initDensityPenalty = initDensityPenalityFactor_;
    npVars.initWireLengthCoef = initWireLengthCoef_;
    npVars.targetOverflow = overflow_;
    npVars.maxNesterovIter = nesterovPlaceMaxIter_;
    npVars.timingDrivenMode = timingDrivenMode_;
    npVars.routabilityDrivenMode = routabilityDrivenMode_;

    std::unique_ptr<NesterovPlace> np(
        new NesterovPlace(npVars, pb_, nb_, log_));
    np_ = std::move(np);
  }
  return true;
}
int Replace::doNesterovPlace(int start_iter) {
  if (!initNesterovPlace()) {
    return 0;
  }
/*
  if (timingDrivenMode_)
    rs_->resizeSlackPreamble();
*/
  return np_->doNesterovPlace(start_iter);
}
void Replace::setInitialPlaceMaxIter(int iter) {
  initialPlaceMaxIter_ = iter;
}

void Replace::setInitialPlaceMinDiffLength(int length) {
  initialPlaceMinDiffLength_ = length;
}

void Replace::setInitialPlaceMaxSolverIter(int iter)
{
  initialPlaceMaxSolverIter_ = iter;
}

void Replace::setInitialPlaceMaxFanout(int fanout)
{
  initialPlaceMaxFanout_ = fanout;
}

void Replace::setInitialPlaceNetWeightScale(float scale)
{
  initialPlaceNetWeightScale_ = scale;
}

void Replace::setNesterovPlaceMaxIter(int iter)
{
  nesterovPlaceMaxIter_ = iter;
  if (np_) {
    np_->setMaxIters(iter);
  }
}

void Replace::setBinGridCntX(int binGridCntX)
{
  binGridCntX_ = binGridCntX;
}

void Replace::setBinGridCntY(int binGridCntY)
{
  binGridCntY_ = binGridCntY;
}

void Replace::setTargetOverflow(float overflow)
{
  overflow_ = overflow;
  if (np_) {
    np_->setTargetOverflow(overflow);
  }
}

void Replace::setTargetDensity(float density)
{
  density_ = density;
}

void Replace::setUniformTargetDensityMode(bool mode)
{
  uniformTargetDensityMode_ = mode;
}

float Replace::getUniformTargetDensity()
{
  initNesterovPlace();
  return nb_->uniformTargetDensity();
}

void Replace::setInitDensityPenalityFactor(float penaltyFactor)
{
  initDensityPenalityFactor_ = penaltyFactor;
}

void Replace::setInitWireLengthCoef(float coef)
{
  initWireLengthCoef_ = coef;
}

void Replace::setMinPhiCoef(float minPhiCoef)
{
  minPhiCoef_ = minPhiCoef;
}

void Replace::setMaxPhiCoef(float maxPhiCoef)
{
  maxPhiCoef_ = maxPhiCoef;
}

void Replace::setReferenceHpwl(float refHpwl)
{
  referenceHpwl_ = refHpwl;
}

void Replace::setSkipIoMode(bool mode)
{
  skipIoMode_ = mode;
}

void Replace::setForceCPU(bool force_cpu)
{
  forceCPU_ = force_cpu;
}

void Replace::setPadLeft(int pad)
{
  padLeft_ = pad;
}

void Replace::setPadRight(int pad)
{
  padRight_ = pad;
}
void Replace::setTimingDrivenMode(bool mode) {
  timingDrivenMode_ = mode;
}

} // gpl