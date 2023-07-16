///////////////////////////////////////////////////////////////////////////////
// Creator: Minjae Kim of CSDL, POSTECH
// Email:   kmj0824@postech.ac.kr
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
#ifndef INC_3D_PLACEMENT_INCLUDE_ALGORITHM_HIERRTLMP_H_
#define INC_3D_PLACEMENT_INCLUDE_ALGORITHM_HIERRTLMP_H_
#include "mpl2/rtl_mp.h"
#include "hier_rtlmp.h"
#include "object.h"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#include "db_sta/MakeDbSta.hh"
#include "sta/Liberty.hh"
#include "par/PartitionMgr.h"
#include "utl/Logger.h"

namespace flow3D {

class HierRTLMPartition : public mpl2::HierRTLMP {
  struct HyperParameters {
    int max_num_macro = 0;
    int min_num_macro = 0;
    int max_num_inst = 0;
    int min_num_inst = 0;
    float tolerance = 0.1;
    int max_num_level = 2;
    float coarsening_ratio = 10.0;
    int num_bundled_ios = 3;
    int large_net_threshold = 50;
    int signature_net_threshold = 50;
    float halo_width = 0.0;
    float fence_lx = 0.0;
    float fence_ly = 0.0;
    float fence_ux = 100000000.0;
    float fence_uy = 100000000.0;
    float area_weight = 0.1;
    float outline_weight = 100.0;
    float wirelength_weight = 100.0;
    float guidance_weight = 10.0;
    float fence_weight = 10.0;
    float boundary_weight = 20.0;
    float notch_weight = 10.0;
    float macro_blockage_weight = 10.0;
    float pin_access_th = 0.00;
    float target_util = 0.25;
    float target_dead_space = 0.05;
    float min_ar = 0.33;
    int snap_layer = -1;
    std::string report_directory = "hier_rtlmp";
  };
 public:
  HierRTLMPartition(sta::dbNetwork *network,
                    odb::dbDatabase *db,
                    sta::dbSta *sta,
                    utl::Logger *logger,
                    par::PartitionMgr *tritonpart) : HierRTLMP(network, db, sta, logger, tritonpart) {
  }
  ~HierRTLMPartition() {
    delete sta_;
    sta::Sta::setSta(nullptr);
  }
  void init();
  void partition();

 private:
  void constructNetwork() {
  }
  void constructSTA() {
  }
  void constructTriton() {
  }

  HyperParameters hyper_parameters_;
};
}
#endif //INC_3D_PLACEMENT_INCLUDE_ALGORITHM_HIERRTLMP_H_
