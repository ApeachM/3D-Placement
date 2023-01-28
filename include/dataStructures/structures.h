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

namespace VLSI_backend {
using namespace odb;
class Instance;
class Net;
class Pin;

/// data storages
struct data_storage {
  std::vector<Instance> instances;
  std::vector<Net> nets;
  std::vector<Pin> pins;
  std::vector<Die> dies;
};

/// data mapping from db to data_storage
struct data_mapping {
  std::unordered_map<dbInst *, Instance *> inst_map;
  std::unordered_map<dbNet *, Net *> net_map;
  /// mapping for terminals on instance (pins on cell)
  std::unordered_map<dbITerm *, Pin *> pin_map_i;
  /// mapping for terminals on blocks (includes fixed pins on die)
  std::unordered_map<dbBTerm *, Pin *> pin_map_b;
};
}
#endif //PLACER_INCLUDE_DATASTRUCTURES_STRUCTURES_H_
