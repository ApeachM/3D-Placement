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
#include "NesterovPlacer.h"
#include "InitialPlacer.h"
#include "Chip.h"

using namespace std;

namespace VLSI_backend {
void Chip::setTargetDensity(vector<double> densities) {
  if (densities.size() != die_pointers_.size())
    assert(0);
  for (int i = 0; i < densities.size(); ++i) {
    die_pointers_.at(i)->setDensity(densities.at(i));
  }
}
void Chip::doInitialPlace() {
  // This function is from below link
  // https://github.com/The-OpenROAD-Project/OpenROAD/blob/977c0794af50e0d3ed993d324b0adead87e32782/src/gpl/src/initialPlace.cpp#L82
  InitialPlacer initial_placer(
      this->instance_pointers_,
      this->net_pointers_,
      this->pin_pointers_,
      this->pad_pointers_,
      this->die_pointers_);

  initial_placer.placeInstancesCenter();
  initial_placer.setPlaceIDs();
  for (int iter = 0; iter < initial_placer.max_iter_; ++iter) {
    initial_placer.updatePinInfo();
    initial_placer.createSparseMatrix();
    pair<float, float> error = initial_placer.cpuSparseSolve();
    float error_max = max(error.first, error.second);
    cout << "[InitialPlace] Iter: " << iter << "\tHPWL: " << getHPWL() << " CG residual: " << error_max << endl;
    if (error_max < 1e-5 && iter >= 5)
      break;
  }
}
void Chip::doNestrovPlace() {
  NesterovPlacer nestrov_placer(
      this->db_database_,
      this->instance_pointers_,
      this->net_pointers_,
      this->pin_pointers_,
      this->pad_pointers_,
      this->die_pointers_.at(0));
  nestrov_placer.initNestrovPlace();
  nestrov_placer.setMaxNesterovIter(30);
  nestrov_placer.doNestrovPlace(0);
}
int Chip::getInstanceNumber() const {
  return instance_number_;
}
void Chip::setInstanceNumber(int instance_number) {
  data_storage_.instances.reserve(instance_number);
  instance_number_ = instance_number;
}
int Chip::getNetNumber() const {
  return net_number_;
}
void Chip::setNetNumber(int net_number) {
  data_storage_.nets.reserve(net_number);
  net_number_ = net_number;
}
dbDatabase *Chip::getDbDatabase() const {
  return db_database_;
}
void Chip::setDbDatabase(dbDatabase *db_database) {
  parser_.db_database_ = db_database;
  db_database_ = db_database;
}

} // VLSI_backend
