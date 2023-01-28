///////////////////////////////////////////////////////////////////////////////
// Creator: Minjae Kim of CSDL, POSTECH
// Email:   kmj0824@postech.ac.kr
// GitHub:  ApeachM
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

#ifndef PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#define PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
#include <utility>
#include <vector>
#include <unordered_map>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include "Parser.h"
#include "Instance.h"
#include "Net.h"
#include "Pin.h"
#include "Die.h"

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;
typedef Eigen::Triplet<float> T;
using Eigen::BiCGSTAB;
using Eigen::IdentityPreconditioner;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;

namespace VLSI_backend {
using namespace odb;

/**\brief
 * This will do 3d placement (D2D placement). \n
 *
 * \usage
 * You should call the methods following below order in normal way.\n
 * Core method is \c do3DPlace(). \n
 * <ol><li> \c parse() \n
 * <li> \c do3DPlace() \n
 * <li> \c write() \n
 *
 * \details
 * In the \c do3DPlace(), placement, partitioning, and synchronized placing will be conducted only one \c odb.\n
 * After all processing, two \c odb(database) will be generated when writing two defs.
 * */
class Chip {
 protected:
  Parser parser_;
  data_storage data_storage_;
  data_mapping data_mapping_;

  std::vector<Instance *> instance_pointers_;
  std::vector<Net *> net_pointers_;
  std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
  std::vector<Pin *> pad_pointers_;
  std::vector<Die *> die_pointers_;

  int num_technologies_ = 0;
  int lib_cell_num_ = 0;
  int util_ = 100;

  /* Placement */
  class InitialPlacer {
   private:
    int max_fan_out_ = 200;
    float net_weight_scale_ = 800.0;
    int min_diff_length_ = 1500;
    int max_solver_iter_ = 100;
    std::vector<Instance *> instance_pointers_;
    std::vector<Net *> net_pointers_;
    std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
    std::vector<Pin *> pad_pointers_;
    std::vector<Die *> die_pointers_;
   public:
    int max_iter_ = 10;  // in OpenROAD, the default value is 20
    InitialPlacer(std::vector<Instance *> instance_pointers,
                  std::vector<Net *> net_pointers,
                  std::vector<Pin *> pin_pointers,
                  std::vector<Pin *> pad_pointers,
                  std::vector<Die *> die_pointers) {
      // TODO: should check move whether method is proper or not.
      instance_pointers_ = std::move(instance_pointers);
      net_pointers_ = std::move(net_pointers);
      pin_pointers_ = std::move(pin_pointers);
      pad_pointers_ = std::move(pad_pointers);
      die_pointers_ = std::move(die_pointers);
    }
    Eigen::VectorXf instLocVecX_, fixedInstForceVecX_;
    Eigen::VectorXf instLocVecY_, fixedInstForceVecY_;
    SMatrix placeInstForceMatrixX_, placeInstForceMatrixY_;
    void placeInstancesCenter() {
      const int center_x = floor(die_pointers_.at(0)->getWidth() / 2);
      const int center_y = floor(die_pointers_.at(0)->getHeight() / 2);

      for (Instance *instance : instance_pointers_) {
        instance->setCoordinate(center_x, center_y);
      }
    }
    void updatePinInfo() {
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
    void setPlaceIDs() {
      for (int i = 0; i < instance_pointers_.size(); ++i) {
        Instance *instance = instance_pointers_.at(i);
        instance->setId(i);
      }
    }
    void createSparseMatrix() {
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
    pair<float, float> cpuSparseSolve() {
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
  };

  void setTargetDensity(vector<double> densities);
  void doInitialPlace() {
    // This function is from below link
    // https://github.com/The-OpenROAD-Project/OpenROAD/blob/977c0794af50e0d3ed993d324b0adead87e32782/src/gpl/src/initialPlace.cpp#L82
    InitialPlacer initial_placer(
        this->instance_pointers_,
        this->net_pointers_,
        this->pin_pointers_,
        this->pad_pointers_,
        this->die_pointers_);

    initial_placer.placeInstancesCenter();
    initial_placer.setPlaceIDs();
    for (int iter = 0; iter < initial_placer.max_iter_; ++iter) {
      initial_placer.updatePinInfo();
      initial_placer.createSparseMatrix();
      pair<float, float> error = initial_placer.cpuSparseSolve();
      float error_max = max(error.first, error.second);
      cout << "[InitialPlace] Iter: " << iter << "HPWL: " << getHPWL() << " CG residual: " << error_max << endl;
      if (error_max < 1e-5 && iter >= 5)
        break;
    }
  }

  /**\brief
   * Placement before partitioning
   * */
  void normalPlacement();

  /*!
 * \brief
 * Do placement 2 Die synchronously.\n
 * This will consider the interaction between two die. \n
 * \details
 * This function is called in `do3DPlace()`.
 * */
  void placement2DieSynchronously();

  /*!
   * \brief
   * Divide a cells into two circuit.
   * Louvain(actually, not louvain but ledien) clustering is implemented by igraph package
   * \todo
   * This code is just temporal code now. Meaningless but only simple partition is implemented.
   * */
  void partition();

  /**\brief
   * get unit of micro
   * \details
   * the coordinate in this circuit is `return value`/1um.
   * \example
   * if the return value is 100, then
   * (20000, 30000) means coordinate (200um, 300um)
   * */
  int getUnitOfMicro() const;
  ulong getHPWL();

  // placer
  void init();

 public:
  Chip() = default;
  ~Chip() = default;
  /*!
   * \brief
   * Core method
   * \details
   * 1. do replace in pseudo die \n
   * 2. partition into two dies
   * 3. do placement synchronously
   * */
  void do3DPlace();

  /*!
   * \brief
   * Read and Write functions
   * */
  void parse(const string &lef_name, const string &def_name);
  void parse_iccad(const string &lef_name, const string &def_name);
  void write(const string &out_file_name);

  // etc
  void dbTutorial() const;
};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
