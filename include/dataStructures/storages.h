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
#ifndef PLACER_INCLUDE_DATASTRUCTURES_STRUCTURES_H_
#define PLACER_INCLUDE_DATASTRUCTURES_STRUCTURES_H_
#include <unordered_map>
#include <vector>
#include "db.h"
#include "Die.h"

namespace flow3D {
using namespace odb;
using namespace std;
class Instance;
class Net;
class Pin;
class HybridBond;

/// data storages
struct data_storage {
  std::vector<Instance> instances;
  std::vector<Net> nets;
  std::vector<Pin> pins;
  std::vector<Die> dies;

  // for hybrid bonds
  std::vector<HybridBond> hybrid_bonds; // set this object as Instance for placing compatibility
};

// ICCAD contest benchmark info
struct LibPinInfo {
  string pin_name;
  int pin_location_x;
  int pin_location_y;
};
struct LibCellInfo {
  string name;
  int width;
  int height;
  int pin_number;
  vector<LibPinInfo> lib_pin_infos;
};
struct TechInfo {
  string name;
  int lib_cell_num;
  vector<LibCellInfo> lib_cell_infos;
};
struct DieInfo {
  string tech_name;
  int lower_left_x;
  int lower_left_y;
  int upper_right_x;
  int upper_right_y;
  int max_util;
  RowInfo row_info;
  TechInfo *tech_info;
};
struct TerminalInfo {
  int size_x;
  int size_y;
  int spacing_size;
};
struct InstanceInfo {
  string inst_name;
  string lib_cell_name;
};
struct ConnectedPinInfo {
  string instance_name;
  string lib_pin_name;
};
struct NetInfo {
  string net_name;
  int pin_num;
  vector<ConnectedPinInfo> connected_pins_infos;
};
struct BenchInformation{
  vector<TechInfo> tech_infos;
  vector<InstanceInfo> instance_infos;
  vector<NetInfo> net_infos;
  vector<DieInfo> die_infos;
  TerminalInfo terminal_info{};
};
enum BENCH_FORMAT {
  ICCAD,
  STANDARD
};
enum BENCH_TYPE{
  ICCAD22,
  ICCAD23,
  ISPD18,
  ETC
};

}
#endif //PLACER_INCLUDE_DATASTRUCTURES_STRUCTURES_H_
