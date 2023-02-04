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
#include <random>
#include <algorithm>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include "Parser.h"
#include "Instance.h"
#include "Net.h"
#include "Pin.h"
#include "Die.h"
#include "fft.h"

#define REPLACE_SQRT2 1.414213562373095048801L

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;
typedef Eigen::Triplet<float> T;
using Eigen::BiCGSTAB;
using Eigen::IdentityPreconditioner;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix;

namespace VLSI_backend {
using namespace odb;
// https://codingforspeed.com/using-faster-exponential-approximation/
static float fastExp(float a);

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
    int max_iter_ = 0;  // in OpenROAD, the default value is 20
    InitialPlacer(std::vector<Instance *> instance_pointers,
                  std::vector<Net *> net_pointers,
                  std::vector<Pin *> pin_pointers,
                  std::vector<Pin *> pad_pointers,
                  std::vector<Die *> die_pointers) {
      instance_pointers_ = instance_pointers;
      net_pointers_ = net_pointers;
      pin_pointers_ = pin_pointers;
      pad_pointers_ = pad_pointers;
      die_pointers_ = die_pointers;
    }
    Eigen::VectorXf instLocVecX_, fixedInstForceVecX_;
    Eigen::VectorXf instLocVecY_, fixedInstForceVecY_;
    SMatrix placeInstForceMatrixX_, placeInstForceMatrixY_;
    void placeInstancesCenter() {
      const int center_x = floor(die_pointers_.at(0)->getWidth() / 2);
      const int center_y = floor(die_pointers_.at(0)->getHeight() / 2);

      for (Instance *instance : instance_pointers_) {
        if (!instance->isLocked())
          instance->setCoordinate(center_x - (instance->getWidth() / 2), center_y - (instance->getHeight() / 2));
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
  class NestrovPlacer {
    class Bin {
     private:
      // index
      int x_;
      int y_;

      // coordinate
      int lx_;
      int ly_;
      int ux_;
      int uy_;

      int64_t nonPlaceArea_;
      int64_t instPlacedArea_;

      int64_t instPlacedAreaUnscaled_;
      int64_t nonPlaceAreaUnscaled_;
      int64_t fillerArea_;

      float density_;
      float targetDensity_;  // will enable bin-wise density screening
      float electroPhi_;
      float electroForceX_;
      float electroForceY_;

     public:
      Bin();
      Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity);

      ~Bin();

      int x() const;
      int y() const;

      int lx() const;
      int ly() const;
      int ux() const;
      int uy() const;
      int cx() const;
      int cy() const;
      int dx() const;
      int dy() const;

      float electroPhi() const;
      float electroForceX() const;
      float electroForceY() const;
      float targetDensity() const;
      float density() const;

      void setDensity(float density);
      void setTargetDensity(float density);
      void setElectroForce(float electroForceX, float electroForceY);
      void setElectroPhi(float phi);

      void setNonPlaceArea(int64_t area);
      void setInstPlacedArea(int64_t area);
      void setFillerArea(int64_t area);

      void setNonPlaceAreaUnscaled(int64_t area);
      void setInstPlacedAreaUnscaled(int64_t area);

      void addNonPlaceArea(int64_t area);
      void addInstPlacedArea(int64_t area);
      void addFillerArea(int64_t area);

      void addNonPlaceAreaUnscaled(int64_t area);
      void addInstPlacedAreaUnscaled(int64_t area);

      const int64_t binArea() const;
      const int64_t nonPlaceArea() const { return nonPlaceArea_; }
      const int64_t instPlacedArea() const { return instPlacedArea_; }
      const int64_t nonPlaceAreaUnscaled() const { return nonPlaceAreaUnscaled_; }
      const int64_t instPlacedAreaUnscaled() const { return instPlacedAreaUnscaled_; }

      const int64_t fillerArea() const { return fillerArea_; }
    };
    class biNormalParameters {
     public:
      float meanX;
      float meanY;
      float sigmaX;
      float sigmaY;
      float lx;
      float ly;
      float ux;
      float uy;
    };

   private:
    odb::dbDatabase *db_database_;
    std::vector<Instance *> instance_pointers_;
    std::vector<Instance *> dummyInsts_;
    std::vector<Instance *> nonPlaceInsts_;
    std::vector<Net *> net_pointers_;
    std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
    std::vector<Pin *> pad_pointers_;
    std::vector<Bin *> bins_;
    Die *die_pointer_ = nullptr;

    // real data storage
    std::vector<Instance> fillers_;
    std::vector<Bin> binStor_;

    // class variables for nestrov place
    int maxNesterovIter = 5000;
    int maxBackTrack = 10;
    float initDensityPenalty = 0.00008;           // INIT_LAMBDA
    float initWireLengthCoef = 0.25;           // base_wcof
    float targetOverflow = 0.1;               // overflow
    float minPhiCoef = 0.95;                   // pcof_min
    float maxPhiCoef = 1.05;                   // pcof_max
    float minPreconditioner = 1.0;            // MIN_PRE
    float initialPrevCoordiUpdateCoef = 100;  // z_ref_alpha
    float referenceHpwl = 446000000;                // refDeltaHpwl
    float routabilityCheckOverflow = 0.20;

    static const int maxRecursionWlCoef{10};
    static const int maxRecursionInitSLPCoef{10};

    int filler_width_{}, filler_height_{};

    int64_t place_instances_area_ = 0;
    int64_t non_place_instances_area_ = 0;

    int64_t std_instances_area_ = 0;
    int64_t macro_instances_area_ = 0;

    int64_t white_space_area_ = 0;
    int64_t movable_area_ = 0;
    int64_t total_filler_area_ = 0;

    float sumPhi_{};
    float target_density_{};
    float uniform_target_density_{};

    bool use_uniform_target_density_ = false;
    bool seed_fix = true;

    int binCntX_ = 0;
    int binCntY_ = 0;
    int binSizeX_ = 0;
    int binSizeY_ = 0;

    int64_t overflowArea_ = 0;
    int64_t overflowAreaUnscaled_ = 0;

    float minWireLengthForceBar = -300;

    gpl::FFT *fft_{};

    // SLP is Step Length Prediction.
    //
    // y_st, y_dst, y_wdst, w_pdst
    std::vector<pair<float, float>> curSLPCoordi_;
    std::vector<pair<float, float>> curSLPWireLengthGrads_;
    std::vector<pair<float, float>> curSLPDensityGrads_;
    std::vector<pair<float, float>> curSLPSumGrads_;

    // y0_st, y0_dst, y0_wdst, y0_pdst
    std::vector<pair<float, float>> nextSLPCoordi_;
    std::vector<pair<float, float>> nextSLPWireLengthGrads_;
    std::vector<pair<float, float>> nextSLPDensityGrads_;
    std::vector<pair<float, float>> nextSLPSumGrads_;

    // z_st, z_dst, z_wdst, z_pdst
    std::vector<pair<float, float>> prevSLPCoordi_;
    std::vector<pair<float, float>> prevSLPWireLengthGrads_;
    std::vector<pair<float, float>> prevSLPDensityGrads_;
    std::vector<pair<float, float>> prevSLPSumGrads_;

    // x_st and x0_st
    std::vector<pair<float, float>> curCoordi_;
    std::vector<pair<float, float>> nextCoordi_;

    // save initial coordinates -- needed for RD
    std::vector<pair<float, float>> initCoordi_;

    // densityPenalty stor
    std::vector<float> densityPenaltyStor_;

    float wireLengthGradSum_{};
    float densityGradSum_{};

    // alpha
    float stepLength_{};

    // opt_phi_cof
    float densityPenalty_{};

    // base_wcof
    float baseWireLengthCoef_{};

    // wlen_cof
    float wireLengthCoefX_{};
    float wireLengthCoefY_{};

    // phi is described in ePlace paper.
    float sumOverflow_{};
    float sumOverflowUnscaled_{};

    // half-parameter-wire-length
    int64_t prevHpwl_{};

    string divergeMsg_;
    float isDiverged_{false};
    float isRoutabilityNeed_{};

    int divergeCode_{};

    int recursionCntWlCoef_{};
    int recursionCntInitSLPCoef_{0};

    bool is_base_initialized = false;

   public:
    NestrovPlacer(odb::dbDatabase *db_database,
                  std::vector<Instance *> instance_pointers,
                  std::vector<Net *> net_pointers,
                  std::vector<Pin *> pin_pointers,
                  std::vector<Pin *> pad_pointers,
                  Die *die_pointer);
    bool initNestrovPlace();
    /*!
     * \brief
     * Do Nestrov Placement. The core function of this class
     * \return
     * the iteration number will be returned by this function
     */

    int doNestrovPlace(int start_iter);;
   private:
    void setInstancesArea();
    void initFillerCells();
    void initBins();
    void updateBinsNonPlaceArea();
    float getDensityCoordiLayoutInsideX(Instance *instance, float cx);
    float getDensityCoordiLayoutInsideY(Instance *instance, float cy);
    void updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates);
    std::pair<int, int> getMinMaxIdxX(Instance *inst) const;
    std::pair<int, int> getMinMaxIdxY(Instance *inst) const;
    // TODO: need to be examined
    int lx() const { return die_pointer_->getLowerLeftX(); }
    int ly() const { return die_pointer_->getLowerLeftY(); }
    int ux() const { return die_pointer_->getUpperRightX(); }
    int uy() const { return die_pointer_->getUpperRightY(); }
    static int fastModulo(int input, const int ceil);
    static int64_t getOverlapArea(const Bin *bin, Instance *inst, int dbu_per_micron);
    static float calculateBiVariateNormalCDF(biNormalParameters i);
    static int64_t getOverlapAreaUnscaled(const Bin *bin, Instance *inst);
    void updateDensitySize();
    void updateDensityForceBin();
    void updateWireLengthForceWA(float wlCoeffX, float wlCoeffY);
    float nesterovInstsArea() const;
    void updateWireLengthCoef(float overflow);
    float getPhiCoef(float scaledDiffHpwl);
    int64_t getHpwl();
    void updateNextIter();
    pair<float, float> getWireLengthPreconditioner(Instance *instance);
    pair<float, float> getDensityPreconditioner(Instance *gCell);
    std::pair<int, int> getDensityMinMaxIdxX(Instance *gcell);
    std::pair<int, int> getDensityMinMaxIdxY(Instance *gcell);
    float getOverlapDensityArea(Bin *bin, Instance *cell);
    pair<float, float> getDensityGradient(Instance *gCell);
    void updateGradients(vector<pair<float, float>> &sumGrads,
                         vector<pair<float, float>> &wireLengthGrads,
                         vector<pair<float, float>> &densityGrads);
    float getStepLength();
    static float getDistance(vector<pair<float, float>> a,
                             vector<pair<float, float>> b);
    pair<float, float> getWireLengthGradientWA(Instance *gCell, float wlCoeffX, float wlCoeffY) const;
    static pair<float, float> getWireLengthGradientPinWA(Pin *gPin, float wlCoeffX, float wlCoeffY) ;
    void updateBinsGCellDensityArea(vector<Instance *> vector_1);
    void setDensityValuesAsDefault();
    void updateDensityCoordiLayoutInside(Instance *gCell);
    void updateInitialPrevSLPCoordi();
    void initSLPStepsVars();
  };

  void setTargetDensity(vector<double> densities);
  void doInitialPlace();
  void doNestrovPlace();

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
