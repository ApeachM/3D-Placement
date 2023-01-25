/*!
 * This code is from https://github.com/The-OpenROAD-Project/OpenROAD
 * */

#ifndef INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_SRC_ALGORITHM_REPLACE_REPLACE_H_
#define INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_SRC_ALGORITHM_REPLACE_REPLACE_H_
#include <memory>
#include "db.h"
#include "initialPlace.h"
//#include "nesterovBase.h"
//#include "nesterovPlace.h"
#include "placerBase.h"
//#include "routeBase.h"
//#include "rsz/Resizer.hh"
//#include "timingBase.h"
//#include "utl/Logger.h"

class PlacerBase;
class NesterovBase;
class InitialPlace;
class NesterovPlace;

namespace gpl {

class Replace {
 public:
  Replace();
  virtual ~Replace();

  void init ();
  void reset();

  void setDb(odb::dbDatabase* odb);

  void doIncrementalPlace();
  void doInitialPlace();

  bool initNesterovPlace();
  int doNesterovPlace(int start_iter = 0);

  // Initial Place param settings
  void setInitialPlaceMaxIter(int iter);
  void setInitialPlaceMinDiffLength(int length);
  void setInitialPlaceMaxSolverIter(int iter);
  void setInitialPlaceMaxFanout(int fanout);
  void setInitialPlaceNetWeightScale(float scale);

  void setNesterovPlaceMaxIter(int iter);

  void setBinGridCntX(int binGridCntX);
  void setBinGridCntY(int binGridCntY);

  void setTargetDensity(float density);
  void setUniformTargetDensityMode(bool mode);
  void setTargetOverflow(float overflow);
  void setInitDensityPenalityFactor(float penaltyFactor);
  void setInitWireLengthCoef(float coef);
  void setMinPhiCoef(float minPhiCoef);
  void setMaxPhiCoef(float maxPhiCoef);

  float getUniformTargetDensity();

  // HPWL: half-parameter wire length.
  void setReferenceHpwl(float deltaHpwl);

  // temp funcs; OpenDB should have these values.
  void setPadLeft(int padding);
  void setPadRight(int padding);

  void setForceCPU(bool force_cpu);
  void setTimingDrivenMode(bool mode);

  void setSkipIoMode(bool mode);


 private:
  odb::dbDatabase * db_;
  utl::Logger* log_;

  std::shared_ptr<PlacerBase> pb_;
  std::shared_ptr<NesterovBase> nb_;

  std::unique_ptr<InitialPlace> ip_;
  std::unique_ptr<NesterovPlace> np_;

  int initialPlaceMaxIter_;
  int initialPlaceMinDiffLength_;
  int initialPlaceMaxSolverIter_;
  int initialPlaceMaxFanout_;
  float initialPlaceNetWeightScale_;
  bool forceCPU_;

  int nesterovPlaceMaxIter_;
  int binGridCntX_;
  int binGridCntY_;
  float overflow_;
  float density_;
  float initDensityPenalityFactor_;
  float initWireLengthCoef_;
  float minPhiCoef_;
  float maxPhiCoef_;
  float referenceHpwl_;

  float routabilityCheckOverflow_;
  float routabilityMaxDensity_;
  float routabilityTargetRcMetric_;
  float routabilityInflationRatioCoef_;
  float routabilityMaxInflationRatio_;

  // routability RC metric coefficients
  float routabilityRcK1_, routabilityRcK2_, routabilityRcK3_, routabilityRcK4_;

  int routabilityMaxBloatIter_;
  int routabilityMaxInflationIter_;

  float timingNetWeightMax_;

  bool timingDrivenMode_;
  bool routabilityDrivenMode_;
  bool uniformTargetDensityMode_;
  bool skipIoMode_;

  std::vector<int> timingNetWeightOverflows_;

  // temp variable; OpenDB should have these values.
  int padLeft_;
  int padRight_;


};

} // gpl

#endif //INC_3D_PLACEMENT_WITH_D2D_VERTICAL_CONNECTIONS_SRC_ALGORITHM_REPLACE_REPLACE_H_
