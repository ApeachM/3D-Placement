#ifndef INC_3D_PLACEMENT_INCLUDE_ALGORITHM_PARTITIONER_H_
#define INC_3D_PLACEMENT_INCLUDE_ALGORITHM_PARTITIONER_H_
#include "Chip.h"
#include "TritonPart.h"

namespace VLSI_backend {
class Chip::Partitioner : par::TritonPart {
 public:
  Partitioner(ord::dbNetwork *network, dbDatabase *db, sta::dbSta *sta, utl::Logger *logger) :
      TritonPart(network, db, sta, logger) {
  }
  virtual ~Partitioner() = default;

  /**
   * Refer to PartitionMgr::tritonPartDesign()
   * */
  void init(const string &design_name);

  void doPartitioning();

  void writeSolution();

 private:
  void readNetList(const string& fixed_file, const string& community_file, const string& group_file){
    if (network_)
      TritonPart::ReadNetlist(fixed_file, community_file, group_file);
    else
      this->ReadNetlist();
  }

  void ReadNetlist();
  string solution_file;
  int dbu;
};
}

#endif //INC_3D_PLACEMENT_INCLUDE_ALGORITHM_PARTITIONER_H_
