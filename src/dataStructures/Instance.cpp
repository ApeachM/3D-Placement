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
#include "Chip.h"
#include "Instance.h"

namespace VLSI_backend {

string Instance::getName() {
  return name_;
}
string Instance::getLibName() {
  return libName_;
}
dbInst *Instance::getDbInst() const {
  if (db_inst_ == nullptr)
    std::cout << "Invalid access to db data from instance" << std::endl;
  return db_inst_;
}
std::vector<Pin *> Instance::getPins() {
  std::vector<Pin *> pins;
  if (db_inst_ != nullptr) {
    // if this is not a filler
    // TODO
    //    can be more simplified?
    for (Pin *pin : connected_pins_) {
      pins.push_back(pin);
    }
  }
  if (is_hybrid_bond_) {
    // TODO
    //    should be need any code?
  }
  return pins;
}
Instance::Instance(odb::dbInst *db_inst) {
  db_database_ = db_inst->getDb();
  db_inst_ = db_inst;
  name_ = db_inst->getName();
  libName_ = db_inst->getMaster()->getName();
  is_macro_ = db_inst_->getMaster()->getType().isBlock();
  position_.first = getCoordinate().first;
  position_.second = getCoordinate().second;
  width_ = db_inst_->getMaster()->getWidth();
  height_ = db_inst_->getMaster()->getHeight();
  dLx_ = position_.first;
  dLy_ = position_.second;
  dUx_ = dLx_ + width_;
  dUy_ = dLy_ + height_;
}
uint Instance::getArea() {
  return this->getWidth() * this->getHeight();
}
void Instance::setCoordinate(int x, int y) {
  position_.first = x;
  position_.second = y;
  if (db_inst_ != nullptr) {
    // if this is normal(not filler) instance,
    db_inst_->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    db_inst_->setLocation(x, y);
  }
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
  if (db_inst_)
    return db_inst_->getMaster()->getWidth();
  else
    // when this instance is filler or hybrid bond
    return width_;
}
uint Instance::getHeight() {
  if (db_inst_)
    return db_inst_->getMaster()->getHeight();
  else
    // when this instance is filler
    return height_;
}
pair<int, int> Instance::getCoordinate() {
  int x = 0, y = 0;
  if (db_inst_)
    db_inst_->getLocation(x, y);
  else {
    x = position_.first;
    y = position_.second;
  }
  return pair<int, int>{x, y};
}
int Instance::getId() const {
  return id_;
}
void Instance::setId(int id) {
  id_ = id;
}
int Instance::dLx() const {
  return dLx_;
}
void Instance::setDLx(int d_lx) {
  dLx_ = d_lx;
}
int Instance::dLy() const {
  return dLy_;
}
void Instance::setDLy(int d_ly) {
  dLy_ = d_ly;
}
int Instance::dUx() const {
  return dUx_;
}
void Instance::setDUx(int d_ux) {
  dUx_ = d_ux;
}
int Instance::dUy() const {
  return dUy_;
}
void Instance::setDUy(int d_uy) {
  dUy_ = d_uy;
}
float Instance::densityScale() const {
  return densityScale_;
}
void Instance::setDensityScale(float density_scale) {
  densityScale_ = density_scale;
}
void Instance::setDensityValueAsDefault() {
  dLx_ = position_.first;
  dLy_ = position_.second;
  dUx_ = position_.first + getWidth();
  dUy_ = position_.second + getHeight();
}
float Instance::getGradientX() const {
  return gradientX_;
}
void Instance::setGradientX(float gradient_x) {
  gradientX_ = gradient_x;
}
float Instance::getGradientY() const {
  return gradientY_;
}
void Instance::setGradientY(float gradient_y) {
  gradientY_ = gradient_y;
}
bool Instance::isMacroInstance() {
  if (!isInstance())
    return false;
  return is_macro_;
}
bool Instance::isStdInstance() {
  if (!isInstance())
    return false;
  return !is_macro_;
}
bool Instance::isFiller() {
  // TODO: make a `is_filler_` variable for this method. We should consider a case for hybrid bond
  return db_inst_ == nullptr;
}
void Instance::setDensityLocation(float dLx, float dLy) {
  dUx_ = dLx + (dUx_ - dLx_);
  dUy_ = dLy + (dUy_ - dLy_);
  dLx_ = dLx;
  dLy_ = dLy;
}
void Instance::setDensityCenterLocation(int dCx, int dCy) {
  const int halfDDx = getDensityDeltaX() / 2;
  const int halfDDy = getDensityDeltaY() / 2;

  dLx_ = dCx - halfDDx;
  dLy_ = dCy - halfDDy;
  dUx_ = dCx + halfDDx;
  dUy_ = dCy + halfDDy;
  for (Pin *pin : getPins()) {
    pin->updateDensityLocation(this);
  }
}
bool Instance::isFixed() {
  // ref: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L139
  // dummy instance is always fixed
  if (isDummy())
    return true;

  switch (db_inst_->getPlacementStatus()) {
    case odb::dbPlacementStatus::NONE:
    case odb::dbPlacementStatus::UNPLACED:
    case odb::dbPlacementStatus::SUGGESTED:
    case odb::dbPlacementStatus::PLACED:return false;
      break;
    case odb::dbPlacementStatus::LOCKED:
    case odb::dbPlacementStatus::FIRM:
    case odb::dbPlacementStatus::COVER:return true;
      break;
  }
  return false;
}
void Instance::setDensitySize(float density_width, float density_height) {
  const uint density_center_x = getDensityCenterX();
  const uint density_center_y = getDensityCenterY();
  dLx_ = static_cast<int>(density_center_x - density_width / 2);
  dLy_ = static_cast<int>(density_center_y - density_height / 2);
  dUx_ = static_cast<int>(density_center_x + density_width / 2);
  dUy_ = static_cast<int>(density_center_y + density_height / 2);
}
int Instance::getDieId() const {
  return die_id_;
}
bool Instance::isHybridBond() const {
  return is_hybrid_bond_;
}
void Instance::setAsHybridBond() {
  is_hybrid_bond_ = true;
  width_ = 0;
  height_ = 0;
  die_id_ = -1;
}
Pin *Instance::getHybridBondPin() const {
  return hybrid_bond_pin_;
}
void Instance::setHybridBondPin(Pin *hybrid_bond_pin) {
  if (!hybrid_bond_pin->isHybridBondPin())
    assert(0); // This pin is not a pin for hybrid bond.
  hybrid_bond_pin_ = hybrid_bond_pin;
}
void Instance::setConnectedPins(vector<Pin *> connected_pins) {
  connected_pins_ = connected_pins;
}
void Instance::setConnectedNets(vector<Net *> connected_nets) {
  connected_nets_ = connected_nets;
}
void Instance::setInstName(const string &name) {
  name_ = name;
}
void Instance::setLibName(const string &lib_name) {
  libName_ = lib_name;
}
dbDatabase *Instance::getDbDatabase() const {
  return db_database_;
}
void Instance::setDbDatabase(dbDatabase *db_database) {
  db_database_ = db_database;
}
void Instance::setLibrary(dbMaster *master) {
  assert(!libName_.empty());
  assert(db_inst_ == nullptr);
  assert(db_database_ == nullptr);
  dbBlock* db_block = db_database_->getChip()->getBlock();
  db_inst_ = dbInst::create(db_block, master, libName_.c_str());
}

} // VLSI_backend
