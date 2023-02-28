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
#ifndef INC_3D_PLACEMENT_INCLUDE_ALGORITHM_INITIALPLACER_H_
#define INC_3D_PLACEMENT_INCLUDE_ALGORITHM_INITIALPLACER_H_
#include "Chip.h"
namespace VLSI_backend {
class Chip::InitialPlacer {
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
}
#endif //INC_3D_PLACEMENT_INCLUDE_ALGORITHM_INITIALPLACER_H_
