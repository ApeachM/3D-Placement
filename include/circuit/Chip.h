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
 *
 * \author
 * Minjae Kim \n
 * GitHub: ApeachM (https://github.com/ApeachM)
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

  /* Placers */
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
                  std::vector<Die *> die_pointers);
    Eigen::VectorXf instLocVecX_, fixedInstForceVecX_;
    Eigen::VectorXf instLocVecY_, fixedInstForceVecY_;
    SMatrix placeInstForceMatrixX_, placeInstForceMatrixY_;
    void placeInstancesCenter();
    void updatePinInfo();
    void setPlaceIDs();
    void createSparseMatrix();
    pair<float, float> cpuSparseSolve();
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

    float curA = 1.0;
    float prevA;

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

    int doNestrovPlace(int start_iter, bool only_one_iter=false);
    void updateDB();

    int getMaxNesterovIter() const;
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

  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void setTargetDensity(vector<double> densities);
  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doInitialPlace();
  /*!
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void doNestrovPlace();

  /**\brief
   * One die (virtual die) placement before partitioning
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void normalPlacement();

  /*!
   * \brief
   * Divide a cells into two circuit.
   * Louvain(actually, not louvain but ledien) clustering is implemented by igraph package
   * \todo
   * This code is just temporal code now. Meaningless but only simple partition is implemented.
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void partition();

  /*!\brief
   * After partitioning, the
   * */
  void generateHybridBonds();

  /*!
   * \brief
   * Do placement 2 Die synchronously.\n
   * This will consider the interaction between two die. \n
   * \details
   * This function is called in `do3DPlace()`.
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
 * */
  void placement2DieSynchronously();

  /**\brief
   * get unit of micro
   * \details
   * the coordinate in this circuit is `return value`/1um.
   * \example
   * if the return value is 100, then
   * (20000, 30000) means coordinate (200um, 300um)
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  int getUnitOfMicro() const;

  /*!
   * \brief
   * get HPWL of total circuit
   * \details
   * This function gets the HPWL by summation of all each nets in `net_pointers`
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  ulong getHPWL();

  // Data initialization
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
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void do3DPlace();

  /*!
   * \brief
   * Read and Write functions
   * \author
   * Minjae Kim \n
   * GitHub: ApeachM (https://github.com/ApeachM)
   * */
  void parse(const string &lef_name, const string &def_name);
  void parse_iccad(const string &lef_name, const string &def_name);
  void write(const string &out_file_name);

  // etc
  void dbTutorial() const;
};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_CIRCUIT_H_
