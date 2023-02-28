#include "InitialPlacer.h"

namespace VLSI_backend {
Chip::InitialPlacer::InitialPlacer(std::vector<Instance *> instance_pointers, std::vector<Net *> net_pointers,
                                   std::vector<Pin *> pin_pointers, std::vector<Pin *> pad_pointers,
                                   std::vector<Die *> die_pointers) {
  instance_pointers_ = instance_pointers;
  net_pointers_ = net_pointers;
  pin_pointers_ = pin_pointers;
  pad_pointers_ = pad_pointers;
  die_pointers_ = die_pointers;
}
void Chip::InitialPlacer::placeInstancesCenter() {
  const int center_x = floor(die_pointers_.at(0)->getWidth() / 2);
  const int center_y = floor(die_pointers_.at(0)->getHeight() / 2);

  for (Instance *instance : instance_pointers_) {
    if (!instance->isLocked())
      instance->setCoordinate(center_x - (instance->getWidth() / 2), center_y - (instance->getHeight() / 2));
  }
}
void Chip::InitialPlacer::updatePinInfo() {
  // set the pin is at the boundary of the bounded box or not

  // reset all MinMax attributes
  for (Pin *pin : pin_pointers_) {
    pin->setMinPinXField(false);
    pin->setMinPinYField(false);
    pin->setMaxPinXField(false);
    pin->setMaxPinYField(false);
  }

  for (Net *net : net_pointers_) {
    Pin *pin_min_x = nullptr, *pin_min_y = nullptr;
    Pin *pin_max_x = nullptr, *pin_max_y = nullptr;

    int lx = INT_MAX, ly = INT_MAX;
    int ux = INT_MIN, uy = INT_MIN;

    // Mark B2B info on Pin structures
    for (Pin *pin : net->getConnectedPins()) {
      if (lx > pin->getCoordinate().first) {
        if (pin_min_x)
          pin_min_x->setMinPinXField(false);
        lx = pin->getCoordinate().first;
        pin_min_x = pin;
        pin_min_x->setMinPinXField(true);
      }

      if (ux < pin->getCoordinate().first) {
        if (pin_max_x)
          pin_max_x->setMaxPinXField(false);
        ux = pin->getCoordinate().first;
        pin_max_x = pin;
        pin_max_x->setMaxPinXField(true);
      }

      if (ly > pin->getCoordinate().second) {
        if (pin_min_y)
          pin_min_y->setMinPinYField(false);
        ly = pin->getCoordinate().second;
        pin_min_y = pin;
        pin_min_y->setMinPinYField(true);
      }

      if (uy < pin->getCoordinate().second) {
        if (pin_max_y)
          pin_max_y->setMaxPinYField(false);
        uy = pin->getCoordinate().second;
        pin_max_y = pin;
        pin_max_y->setMaxPinYField(true);
      }
    }
  }
}
void Chip::InitialPlacer::createSparseMatrix() {
  // This function is from below link
  // https://github.com/The-OpenROAD-Project/OpenROAD/blob/977c0794af50e0d3ed993d324b0adead87e32782/src/gpl/src/initialPlace.cpp#L234
  const int placeCnt = instance_pointers_.size();
  instLocVecX_.resize(placeCnt);
  fixedInstForceVecX_.resize(placeCnt);
  instLocVecY_.resize(placeCnt);
  fixedInstForceVecY_.resize(placeCnt);

  placeInstForceMatrixX_.resize(placeCnt, placeCnt);
  placeInstForceMatrixY_.resize(placeCnt, placeCnt);

  //
  // list_x and list_y is a temporary vector that have tuples, (idx1, idx2, val)
  //
  // list_x finally becomes placeInstForceMatrixX_
  // list_y finally becomes placeInstForceMatrixY_
  //
  // The triplet vector is recommended usages
  // to fill in SparseMatrix from Eigen docs.
  //
  vector<T> list_x, list_y;
  list_x.reserve(1000000);
  list_y.reserve(1000000);

  // initialize vector
  for (Instance *instance : instance_pointers_) {
    int idx = instance->getId();
    instLocVecX_(idx) = instance->getCenterX();
    instLocVecY_(idx) = instance->getCenterY();

    fixedInstForceVecX_(idx) = 0;
    fixedInstForceVecY_(idx) = 0;
  }

  // for each net
  for (Net *net : net_pointers_) {
    // skip for small nets.
    if (net->getConnectedPins().size() <= 1)
      continue;

    // escape long time cals on huge fanout.
    if (net->getConnectedPins().size() >= max_fan_out_)
      continue;

    float net_weight = net_weight_scale_ / static_cast<float>(net->getConnectedPins().size() - 1);

    vector<Pin *> pins = net->getConnectedPins();
    for (int pin_idx1 = 1; pin_idx1 < pins.size(); ++pin_idx1) {
      Pin *pin1 = pins.at(pin_idx1);
      for (int pin_idx2 = 0; pin_idx2 < pin_idx1; ++pin_idx2) {
        Pin *pin2 = pins.at(pin_idx2);

        // no need to fill in when instance is same
        if (pin1->getInstance() == pin2->getInstance())
          continue;

        // B2B modeling on min_x/max_x pins.
        if (pin1->isMinPinX() || pin1->isMaxPinX() || pin2->isMinPinX() || pin2->isMaxPinX()) {
          int diff_x = abs(pin1->getCoordinate().first - pin2->getCoordinate().first);
          float weight_x = 0;
          if (diff_x > min_diff_length_)
            weight_x = net_weight / diff_x;
          else
            weight_x = net_weight / min_diff_length_;

          // both pin came from instance
          if (pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            const int inst2 = pin2->getInstance()->getId();

            list_x.push_back(T(inst1, inst2, weight_x));
            list_x.push_back(T(inst2, inst1, weight_x));
            list_x.push_back(T(inst1, inst2, -weight_x));
            list_x.push_back(T(inst2, inst1, -weight_x));

            fixedInstForceVecX_(inst1) +=
                -weight_x * static_cast<float>
                (
                    (pin1->getCoordinate().first - pin1->getInstance()->getCenterX())
                        - (pin2->getCoordinate().first - pin2->getInstance()->getCenterX())
                );
            fixedInstForceVecX_(inst2) +=
                -weight_x * static_cast<float>
                (
                    (pin2->getCoordinate().first - pin2->getInstance()->getCenterX())
                        - (pin1->getCoordinate().first - pin1->getInstance()->getCenterX())
                );
          }
            // pin1 from IO port / pin2 from Instance
          else if (!pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst2 = pin2->getInstance()->getId();
            list_x.push_back(T(inst2, inst2, weight_x));
            fixedInstForceVecX_(inst2) += weight_x * static_cast<float>
            (pin1->getCoordinate().first - (pin2->getCoordinate().first - pin2->getInstance()->getCenterX()));
          }
            // pin1 from Instance / pin2 from IO port
          else if (pin1->isInstancePin() && !pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            list_x.push_back(T(inst1, inst1, weight_x));
            fixedInstForceVecX_(inst1) += weight_x * static_cast<float>
            (pin2->getCoordinate().first - (pin1->getCoordinate().first - pin1->getInstance()->getCenterX()));
          }
        }

        // B2B modeling on min_y/max_y pins.
        if (pin1->isMinPinY() || pin1->isMaxPinY() || pin2->isMinPinY() || pin2->isMaxPinY()) {
          int diff_y = abs(pin1->getCoordinate().second - pin2->getCoordinate().second);
          float weight_y = 0;
          if (diff_y > min_diff_length_) {
            weight_y = net_weight / static_cast<float>(diff_y);
          } else {
            weight_y = net_weight / min_diff_length_;
          }

          // both pin came from instance
          if (pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            const int inst2 = pin2->getInstance()->getId();
            list_y.push_back(T(inst1, inst1, weight_y));
            list_y.push_back(T(inst2, inst2, weight_y));
            list_y.push_back(T(inst1, inst2, -weight_y));
            list_y.push_back(T(inst2, inst1, -weight_y));

            fixedInstForceVecY_(inst1) += -weight_y * static_cast<float>
            (
                (pin1->getCoordinate().second - pin1->getInstance()->getCenterY())
                    - (pin2->getCoordinate().second - pin2->getInstance()->getCenterY())
            );
            fixedInstForceVecY_(inst2) += -weight_y * static_cast<float>
            (
                (pin2->getCoordinate().second - pin2->getInstance()->getCenterY())
                    - (pin1->getCoordinate().second - pin1->getInstance()->getCenterY())
            );
          }
            // pin1 from IO port // pin2 from Instance
          else if (!pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst2 = pin2->getInstance()->getId();
            list_y.push_back(T(inst2, inst2, weight_y));
            fixedInstForceVecY_(inst2) += weight_y * static_cast<float>(
                (pin1->getCoordinate().second - (pin2->getCoordinate().second - pin2->getInstance()->getCenterY()))
            );
          }
            // pin1 from Instance / pin2 from IO port
          else if (pin1->isInstancePin() && !pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            list_y.push_back(T(inst1, inst1, weight_y));
            fixedInstForceVecY_(inst1) += weight_y * static_cast<float>(
                (pin2->getCoordinate().second - (pin1->getCoordinate().second - pin1->getInstance()->getCenterY()))
            );

          }
        }
      }
    }
  }
  placeInstForceMatrixX_.setFromTriplets(list_x.begin(), list_x.end());
  placeInstForceMatrixY_.setFromTriplets(list_y.begin(), list_y.end());
}
void Chip::InitialPlacer::setPlaceIDs() {
  // reset ExtId for all instances
  for (Instance *instance : instance_pointers_) {
    instance->setId(INT_MAX);
  }
  // set index only with place-able instances
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    Instance *instance = instance_pointers_.at(i);
    if (!instance->isFixed())
      instance->setId(i);
  }
}
pair<float, float> Chip::InitialPlacer::cpuSparseSolve() {
  pair<float, float> error;
  BiCGSTAB<SMatrix, IdentityPreconditioner> solver;
  solver.setMaxIterations(max_solver_iter_);
  // for x
  solver.compute(placeInstForceMatrixX_);
  instLocVecX_ = solver.solveWithGuess(fixedInstForceVecX_, instLocVecX_);
  error.first = solver.error();
  // for y
  solver.compute(placeInstForceMatrixY_);
  instLocVecY_ = solver.solveWithGuess(fixedInstForceVecY_, instLocVecY_);
  error.second = solver.error();

  return error;
}
}