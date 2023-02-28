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
  class NesterovPlacer {
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

    float sum_phi_{};
    float target_density_{};
    float uniform_target_density_{};

    bool use_uniform_target_density_ = false;
    bool seed_fix = true;

    int bin_cnt_x_ = 0;
    int bin_cnt_y_ = 0;
    int bin_size_x_ = 0;
    int bin_size_y_ = 0;

    int64_t overflow_area_ = 0;
    int64_t overflow_area_unscaled_ = 0;

    float min_wire_length_force_bar_ = -300;

    gpl::FFT *fft_{};

    // SLP is Step Length Prediction.
    //
    // y_st, y_dst, y_wdst, w_pdst
    std::vector<pair<float, float>> cur_SLP_coordinates_;
    std::vector<pair<float, float>> cur_SLP_wire_length_grads_;
    std::vector<pair<float, float>> cur_SLP_density_grads_;
    std::vector<pair<float, float>> cur_SLP_sum_grads_;

    // y0_st, y0_dst, y0_wdst, y0_pdst
    std::vector<pair<float, float>> next_SLP_coordinates_;
    std::vector<pair<float, float>> next_SLP_wire_length_grads_;
    std::vector<pair<float, float>> next_SLP_density_grads_;
    std::vector<pair<float, float>> next_SLP_sum_grads_;

    // z_st, z_dst, z_wdst, z_pdst
    std::vector<pair<float, float>> prev_SLP_coordinates_;
    std::vector<pair<float, float>> prev_SLP_wire_length_grads_;
    std::vector<pair<float, float>> prev_SLP_density_grads_;
    std::vector<pair<float, float>> prev_SLP_sum_grads_;

    // x_st and x0_st
    std::vector<pair<float, float>> cur_coordinates_;
    std::vector<pair<float, float>> next_coordinates_;

    // save initial coordinates -- needed for RD
    std::vector<pair<float, float>> init_coordinates_;

    // densityPenalty stor
    std::vector<float> density_penalty_storage_;

    float wire_length_grad_sum_{};
    float density_grad_sum_{};

    // alpha
    float step_length_{};

    // opt_phi_cof
    float density_penalty_{};

    // base_wcof
    float base_wire_length_coefficient_{};

    // wlen_cof
    float wire_length_coefficient_x_{};
    float wire_length_coefficient_y_{};

    // phi is described in ePlace paper.
    float sum_overflow_{};
    float sum_overflow_unscaled_{};

    // half-parameter-wire-length
    int64_t prev_hpwl_{};

    string diverge_msg_;
    float is_diverged_{false};
    float is_routability_need_{};

    int diverge_code_{};

    int recursion_cnt_wl_coef_{};
    int recursion_cnt_init_slp_coef_{0};

    bool is_base_initialized_ = false;

   public:
    NesterovPlacer(odb::dbDatabase *db_database,
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
    /*!
     * \name
     * setInstancesArea
     *
     * \brief
     * Set area counting variables by iterating the all instances
     *
     * \details
     * This method highly refers to https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L798 \n
     * , which is name of `void PlacerBase::init()`.\n\n
     *
     * The set variables are next: `place_instances_area_`, `macro_instances_area_`, and `std_instances_area_`. \n
     * This code should be called at initialization of Nestrov.
     * */
    void setInstancesArea();

    /*!
     * \name
     * initFillerCells
     *
     * \brief
     * This method initialize the fillers.
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L1184 \n
     * , which name is `void NesterovBase::initFillerGCells()`. \n\n
     * This method set the filler width and heights, and the number of the fillers.
     * The number of the fillers considers the area of die, target density, and the area of filler, etc. \n
     * After setting them, this function sets the filler coordinates randomly.
     *
     * \pre
     * You should call `setInstancesArea()` method before calling this code.
     * */
    void initFillerCells();

    /*!
     * \name
     * initBins
     *
     * \brief
     * This methods initializes the bins.
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L723 \n
     * , which name is `void BinGrid::initBins()`. \n\n
     * The \c Bin is for defining the gradient vectors.
     * Using bin notion, the complexity is reasonably highly deduced keeping the result reasonable.\n
     * In this function, the number of bins and the size of bins is determined. \n
     * And in the `updateBinsNonPlaceArea()` function, get the overlapped area with each bin and each instances. \n
     * How many cells is overlapped with the bin determine the charge amounts in the bin.
     *
     * \pre
     * You should call `setInstancesArea()` and `initFillerCells()` methods before calling this code.
     * */
    void initBins();

    /*!
     * \name
     * updateBinsNonPlaceArea
     *
     * \brief
     * Set the non place area for each bin
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L806 \n
     * , which name is `void BinGrid::updateBinsNonPlaceArea()`. \n\n
     * This function should be called only in `initBins()` function. \n
     * Here, the function sets the value for how many cells is overlapped with the each bins.
     *
     * \pre
     * You should call `setInstancesArea()` and `initFillerCells()` methods before calling this code,
     * because this function should be only in `initBins()`.
     * */
    void updateBinsNonPlaceArea();

    /*!
     * \name
     * getDensityCoordiLayoutInsideX
     *
     * \brief
     * If the cell is out of the die,
     * then return the value that the center of the cell coordinate which make the cell be in the die.\n
     * This function is only for x coordinate.
     *
     * \param
     * Instance, float \n
     * float is for the target coordinated that will be adjusted. \n
     * Instance object will be considered when the target coordinate is adjusted.
     *
     * \return
     * float \n
     * adjusted target coordinate
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/e0983e4988d09bcffe31590ae3d921489159fd10/src/gpl/src/nesterovBase.cpp#L1586 \n
     * , which name is `getDensityCoordiLayoutInsideX` . \n\n
     * The `cx`  is considered as target coordinate. This means the center location of the cell. \n
     * If `cx` makes the cell out of the die, then adjust the target coordinate (`adjVal`) to be in the die for the cell.
     * The adjusted target coordinate will be returned as return value. \n
     * This function only for x coordinate.
     *
     * */
    float getDensityCoordiLayoutInsideX(Instance *instance, float cx);

    /*!
     * \name
     * getDensityCoordiLayoutInsideY
     *
     * \brief
     * If the cell is out of the die,
     * then return the value that the center of the cell coordinate which make the cell be in the die.\n
     * This function is only for y coordinate.
     *
     * \param
     * instance This instance object will be considered when the target coordinate is adjusted.
     *
     * \param
     * cy This float variable is for the target coordinated that will be adjusted. \n
     *
     * \return
     * adjusted target coordinate for y axis
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/e0983e4988d09bcffe31590ae3d921489159fd10/src/gpl/src/nesterovBase.cpp#L1586 \n
     * , which name is `getDensityCoordiLayoutInsideX` . \n\n
     * The `cx`  is considered as target coordinate. This means the center location of the cell. \n
     * If `cx` makes the cell out of the die, then adjust the target coordinate (`adjVal`) to be in the die for the cell.
     * The adjusted target coordinate will be returned as return value. \n
     * This function only for y coordinate.
     *
     * */
    float getDensityCoordiLayoutInsideY(Instance *instance, float cy);

    /*!
     * \name
     * updateGCellDensityCenterLocation
     *
     * \brief
     *
     * \param
     * coordinates This coordinates will apply on each cell coordinate for electrical density
     *
     * \details
     * This method highly refers to the code in below link. \n
     * https://github.com/The-OpenROAD-Project/OpenROAD/blob/e0983e4988d09bcffe31590ae3d921489159fd10/src/gpl/src/nesterovBase.cpp#L1352 \n
     * , which name is `updateGCellDensityCenterLocation`. \n\n
     * This function gets the coordinates as the parameter. \n
     * And, the process is divided as two step.\n\n
     * 1. Update the coordinates variables for density things of cells.
     *    The input parameter (the set of coordinates) will be inserted into cell coordinate variables.
     *    \n
     * 2. Call `updateBinsGCellDensityArea` function.
     * This function will considers the cell coordinates variables for density things,
     * and update the bin variables
     * (`instPlacedArea_`, `instPlacedAreaUnscaled_`, `instPlacedAreaUnscaled_`, `nonPlaceAreaUnscaled_`, and `fillerArea_`).
     * */
    void updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates);

    /*!
     *
     * */
    std::pair<int, int> getMinMaxIdxX(Instance *inst) const;
    std::pair<int, int> getMinMaxIdxY(Instance *inst) const;
    // TODO: need to be examined
    /*!
     * \return
     * x coordinate of lower left for die
     * */
    int lx() const { return die_pointer_->getLowerLeftX(); }

    /*!
     * \return
     * y coordinate of lower left for die
     * */
    int ly() const { return die_pointer_->getLowerLeftY(); }
    /*!
     * \return
     * x coordinate of upper right for die
     * */
    int ux() const { return die_pointer_->getUpperRightX(); }
    /*!
     * \return
     * y coordinate of upper right die
     * */
    int uy() const { return die_pointer_->getUpperRightY(); }

    static int fastModulo(int input, const int ceil);
    /*!
     * \return
     * the overlapped area between bin and instance
     * */
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
  class Test;
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
