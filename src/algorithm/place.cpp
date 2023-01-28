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

#include "Chip.h"
#include <random>
namespace VLSI_backend {
void Chip::normalPlacement() {
  /* top die util setting manually in code level */
  vector<double> densities;
  densities.push_back(1.0);
  densities.push_back(1.0);
  setTargetDensity(densities);

  doInitialPlace();


}

void Chip::partition() {
  /* Temporal code */
  int cell_num = static_cast<int>(instance_pointers_.size());
  for (int i = 0; i < floor(cell_num / 2); ++i) {
    Instance* instance = instance_pointers_.at(i);
    instance->assignDie(1);
  }
  for (int i = floor(cell_num/2); i < cell_num; ++i) {
    Instance* instance = instance_pointers_.at(i);
    instance->assignDie(2);
  }
}

void Chip::placement2DieSynchronously() {

}

void Chip::do3DPlace() {
  // 1. do3DPlace the cells in the pseudo die
  this->normalPlacement();

  // 2. partition
  this->partition();

  // 3. placement synchronously
  this->placement2DieSynchronously();
}
}
