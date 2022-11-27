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
#include "Instance.h"

namespace Circuit {

string Instance::getName() {
  return name_;
}
string Instance::getLibName() {
  return libName_;
}
void Instance::setDataStorage(data_storage *data_storage) {
  data_storage_ = data_storage;
}
void Instance::setDataMapping(data_mapping *data_mapping) {
  data_mapping_ = data_mapping;
}
dbInst *Instance::getDbInst() const {
  if(db_inst_ == nullptr)
    std::cout << "Invalid access to db data from instance" << std::endl;
  return db_inst_;
}
std::vector<Pin *> Instance::getPins() {
  std::vector<Pin *> pins;
  for (dbITerm *db_i_term : db_inst_->getITerms()) {
    Pin *pin = data_mapping_->pin_map_i[db_i_term];
    pins.push_back(pin);
  }
  return pins;
}
Instance::Instance(odb::dbInst *db_inst) {
  db_database_ = db_inst->getDb();
  db_inst_ = db_inst;
  name_ = db_inst->getName();
  libName_ = db_inst->getMaster()->getName();
}
Instance::Instance(odb::dbInst *db_inst, data_storage *data_storage, data_mapping *data_mapping) {
  db_database_ = db_inst->getDb();
  data_storage_ = data_storage;
  data_mapping_ = data_mapping;

  db_inst_ = db_inst;
  name_ = db_inst->getName();
  libName_ = db_inst->getMaster()->getName();

}
uint Instance::getArea() {
  return this->getWidth() * this->getHeight();
}
void Instance::setCoordinate(int x, int y) {
  position_.first = x;
  position_.second = y;
  db_inst_->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  db_inst_->setLocation(x, y);
}
bool Instance::isPlaced() {
  if (db_inst_->getPlacementStatus() == odb::dbPlacementStatus::PLACED) {
    return true;
  } else if (db_inst_->getPlacementStatus() == odb::dbPlacementStatus::NONE) {
    return false;
  } else if (db_inst_->getPlacementStatus() == odb::dbPlacementStatus::UNPLACED) {
    return false;
  } else {
    return true;
  }
}
uint Instance::getWidth() {
  return db_inst_->getMaster()->getWidth();
}
uint Instance::getHeight() {
  return db_inst_->getMaster()->getHeight();
}
pair<int, int> Instance::getCoordinate() {
  int x = 0, y = 0;
  db_inst_->getLocation(x, y);
  return pair<int, int>{x, y};
}

} // Circuit
