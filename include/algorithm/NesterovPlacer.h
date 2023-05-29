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
#ifndef INC_3D_PLACEMENT_INCLUDE_ALGORITHM_NESTEROVPLACER_H_
#define INC_3D_PLACEMENT_INCLUDE_ALGORITHM_NESTEROVPLACER_H_
#include "Chip.h"
#include "CImg.h"

namespace VLSI_backend {
class Chip::NesterovPlacer {
  class Bin;
  class biNormalParameters;
  class Drawer;
 public:
  NesterovPlacer(odb::dbDatabase *db_database,
                 std::vector<Instance *> instance_pointers,
                 std::vector<Net *> net_pointers,
                 std::vector<Pin *> pin_pointers,
                 std::vector<Pin *> pad_pointers,
                 Die *die_pointer);
  bool initNesterovPlace(bool is_pseudo_die = true);
  /*!
   * \brief
   * Do Nesterov Placement. The core function of this class
   * \return
   * the iteration number will be returned by this function
   */
  int doNesterovPlace(int start_iter = 0, bool only_one_iter = false);
  void updateDB();

  int getMaxNesterovIter() const { return max_nesterov_iter_; }
  void setMaxNesterovIter(int max_nesterov_iter) { max_nesterov_iter_ = max_nesterov_iter; }
 private:
  /*!
   * \name
   * setInstancesArea
   *
   * \brief
   * Set area counting variables by iterating the all instances
   *
   * \details
   * The set variables are next: `place_instances_area_`, `macro_instances_area_`, and `std_instances_area_`. \n
   * This code should be called at initialization of Nesterov.\n\n
   * This method highly refers to https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L798 \n
   * , which is name of `void PlacerBase::init()`.
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
   * This method set the filler width and heights, and the number of the fillers.
   * The number of the fillers considers the area of die, target density, and the area of filler, etc. \n
   * After setting them, this function sets the filler coordinates randomly.\n\n
   * This method highly refers to the code in below link. \n
   * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L1184 \n
   * , which name is `void NesterovBase::initFillerGCells()`.
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
   * The \c Bin is for defining the gradient vectors.
   * Using bin notion, the complexity is reasonably highly deduced keeping the result reasonable.\n
   * In this function, the number of bins and the size of bins is determined. \n
   * And in the `updateBinsNonPlaceArea()` function, get the overlapped area with each bin and each instances. \n
   * How many cells is overlapped with the bin determine the charge amounts in the bin.\n\n
   * This method highly refers to the code in below link. \n
   * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L723 \n
   * , which name is `void BinGrid::initBins()`.
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
   * This function should be called only in `initBins()` function. \n
   * Here, the function sets the value for how many cells is overlapped with the each bins.\n\n
   * This method highly refers to the code in below link. \n
   * https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L806 \n
   * , which name is `void BinGrid::updateBinsNonPlaceArea()`.
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
   * The `cx`  is considered as target coordinate. This means the center location of the cell. \n
   * If `cx` makes the cell out of the die, then adjust the target coordinate (`adjVal`) to be in the die for the cell.
   * The adjusted target coordinate will be returned as return value. \n
   * This function only for x coordinate. \n\n
   * This method highly refers to the code in below link. \n
   * https://github.com/The-OpenROAD-Project/OpenROAD/blob/e0983e4988d09bcffe31590ae3d921489159fd10/src/gpl/src/nesterovBase.cpp#L1586 \n
   * , which name is `getDensityCoordiLayoutInsideX` .
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
   * 2. Call `updateBinsCellDensityArea` function.
   * This function will considers the cell coordinates variables for density things,
   * and update the bin variables
   * (`instPlacedArea_`, `instPlacedAreaUnscaled_`, `instPlacedAreaUnscaled_`, `nonPlaceAreaUnscaled_`, and `fillerArea_`).
   * */
  void updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates);

  /*!
   * \name
   * getMinMaxIdxX
   * \details
   * */
  std::pair<int, int> getMinMaxIdxX(Instance *inst) const;
  /*!
   * \name
   * getMinMaxIdxY
   * \details
   * */
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
   * \name
   * getOverlapArea
   * \return
   * the overlapped area between bin and instance
   * */
  static int64_t getOverlapArea(const Bin *bin, Instance *inst, int dbu_per_micron);
  static int64_t getOverlapAreaUnscaled(const Bin *bin, Instance *inst);
  /*!
   * \name
   * calculateBiVariateNormalCDF
   * \brief
   *
   * \details
   *   A function that does 2D integration to the density function of a
   *   bivariate normal distribution with 0 correlation.
   *   Essentially, the function being integrated is the product
   *   of 2 1D probability density functions (for x and y). The means and standard
   *   deviation of the probablity density functions are parametarized. In this
   *   function, I am using the closed-form solution of the integration. The limits
   *   of integration are lx->ux and ly->uy For reference: the equation that is
   *   being integrated is:
   *   (1/(2*pi*sigmaX*sigmaY))*e^(-(y-meanY)^2/(2*sigmaY*sigmaY))*e^(-(x-meanX)^2/(2*sigmaX*sigmaX))
   * \pre
   * This is called in only getOverlapArea function
   * */
  static float calculateBiVariateNormalCDF(biNormalParameters i);
  /*!
   * \name
   * updateDensitySize
   *
   * */
  void updateDensitySize();
  void updateDensityForceBin();
  void updateWireLengthForceWA(double wlCoeffX, double wlCoeffY);
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
  static pair<float, float> getWireLengthGradientPinWA(Pin *gPin, float wlCoeffX, float wlCoeffY);
  void updateBinsCellDensityArea(vector<Instance *> cells);
  void setDensityValuesAsDefault();
  void updateDensityCoordiLayoutInside(Instance *gCell);
  void updateInitialPrevSLPCoordi();
  void initSLPStepsVars();
 public:
  void setDebugMode(bool debug_mode);

 private:
  void handleDiverge(const vector<pair<float, float>> &snapshotCoordi,
                     const vector<pair<float, float>> &snapshotSLPCoordi,
                     const vector<pair<float, float>> &snapshotSLPSumGrads,
                     float snapshotA,
                     float snapshotDensityPenalty,
                     float snapshotStepLength,
                     float snapshotWlCoefX,
                     float snapshotWlCoefY,
                     bool &isDivergeTriedRevert);
  bool stepLengthDivergeCheck();
  void printStateNesterov(int iter) const;
  bool finishCheck() const;
  void drawDie(const string &filename);

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

  // class variables for nesterov place
  // int max_nesterov_iter_ = 5000;
  int max_nesterov_iter_ = 200;
  int max_back_track_ = 10;
  float initDensityPenalty = 0.00008;         // INIT_LAMBDA
  float initWireLengthCoef = 0.25;            // base_wcof
  float targetOverflow = 0.1;                 // overflow
  float minPhiCoef = 0.95;                    // pcof_min
  float max_phi_coef_ = 1.05;                    // pcof_max
  float minPreconditioner = 1.0;              // MIN_PRE
  float initialPrevCoordiUpdateCoef = 100;    // z_ref_alpha
  float referenceHpwl = 446000000;            // refDeltaHpwl
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
  std::vector<pair<float, float>> cur_slp_coordinates_;
  std::vector<pair<float, float>> cur_slp_wire_length_grads_;
  std::vector<pair<float, float>> cur_slp_density_grads_;
  std::vector<pair<float, float>> cur_slp_sum_grads_;

  // y0_st, y0_dst, y0_wdst, y0_pdst
  std::vector<pair<float, float>> next_slp_coordinates_;
  std::vector<pair<float, float>> next_slp_wire_length_grads_;
  std::vector<pair<float, float>> next_slp_density_grads_;
  std::vector<pair<float, float>> next_slp_sum_grads_;

  // z_st, z_dst, z_wdst, z_pdst
  std::vector<pair<float, float>> prev_slp_coordinates_;
  std::vector<pair<float, float>> prev_slp_wire_length_grads_;
  std::vector<pair<float, float>> prev_slp_density_grads_;
  std::vector<pair<float, float>> prev_slp_sum_grads_;

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
  bool debug_mode_ = false;
};

class Chip::NesterovPlacer::Bin {
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
class Chip::NesterovPlacer::biNormalParameters {
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
class Chip::NesterovPlacer::Drawer {
  using Image = cimg_library::CImg<unsigned char>;
 public:
  explicit Drawer(uint width = 0, uint height = 0);
  virtual ~Drawer();
  void drawCell(int ll_x, int ll_y, int ur_x, int ur_y);
  void drawFiller(int ll_x, int ll_y, int ur_x, int ur_y);
  void setCellColor(const unsigned char *cell_color);
  void setFillerColor(const unsigned char *filler_color);
  void saveImg(const string &file_name);

 private:
  uint width_;
  uint height_;
  Image *image_;
  const unsigned char *cell_color_ = Color::BLACK;
  const unsigned char *filler_color_ = Color::RED;
  string file_path_ = "../output/images/";
};
}

#endif //INC_3D_PLACEMENT_INCLUDE_ALGORITHM_NESTEROVPLACER_H_
