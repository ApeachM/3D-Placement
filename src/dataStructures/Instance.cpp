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
  d_lx_ = position_.first;
  d_ly_ = position_.second;
  d_ux_ = d_lx_ + width_;
  d_uy_ = d_ly_ + height_;
}
Instance::Instance(odb::dbInst *db_inst, int id) {
  db_database_ = db_inst->getDb();
  db_inst_ = db_inst;
  name_ = db_inst->getName();
  libName_ = db_inst->getMaster()->getName();
  is_macro_ = db_inst_->getMaster()->getType().isBlock();
  position_.first = getCoordinate().first;
  position_.second = getCoordinate().second;
  width_ = db_inst_->getMaster()->getWidth();
  height_ = db_inst_->getMaster()->getHeight();
  d_lx_ = position_.first;
  d_ly_ = position_.second;
  d_ux_ = d_lx_ + width_;
  d_uy_ = d_ly_ + height_;
  setId(id);
}
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
    pins.push_back(hybrid_bond_pin_);
  }
  return pins;
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
  } else if (is_hybrid_bond_) {
    // if this is a hybrid bond, the pin coordinate is not automatically adjusted when instance is moved.
    // so, we should adjust the pin coordinate manually.
    hybrid_bond_pin_->setHybridBondCoordinate(x, y);
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
  return d_lx_;
}
void Instance::setDLx(int d_lx) {
  d_lx_ = d_lx;
}
int Instance::dLy() const {
  return d_ly_;
}
void Instance::setDLy(int d_ly) {
  d_ly_ = d_ly;
}
int Instance::dUx() const {
  return d_ux_;
}
void Instance::setDUx(int d_ux) {
  d_ux_ = d_ux;
}
int Instance::dUy() const {
  return d_uy_;
}
void Instance::setDUy(int d_uy) {
  d_uy_ = d_uy;
}
float Instance::densityScale() const {
  return density_scale_;
}
void Instance::setDensityScale(float density_scale) {
  density_scale_ = density_scale;
}
void Instance::setDensityValueAsDefault() {
  d_lx_ = position_.first;
  d_ly_ = position_.second;
  d_ux_ = position_.first + static_cast<int>(getWidth());
  d_uy_ = position_.second + static_cast<int>(getHeight());
}
float Instance::getGradientX() const {
  return gradient_x_;
}
void Instance::setGradientX(float gradient_x) {
  gradient_x_ = gradient_x;
}
float Instance::getGradientY() const {
  return gradient_y_;
}
void Instance::setGradientY(float gradient_y) {
  gradient_y_ = gradient_y;
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
  d_ux_ = dLx + (d_ux_ - d_lx_);
  d_uy_ = dLy + (d_uy_ - d_ly_);
  d_lx_ = dLx;
  d_ly_ = dLy;
}
void Instance::setDensityCenterLocation(int d_cx, int d_cy) {
  const int half_d_dx = getDensityDeltaX() / 2;
  const int half_d_dy = getDensityDeltaY() / 2;

  d_lx_ = d_cx - half_d_dx;
  d_ly_ = d_cy - half_d_dy;
  d_ux_ = d_cx + half_d_dx;
  d_uy_ = d_cy + half_d_dy;
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
  d_lx_ = static_cast<int>(density_center_x - density_width / 2);
  d_ly_ = static_cast<int>(density_center_y - density_height / 2);
  d_ux_ = static_cast<int>(density_center_x + density_width / 2);
  d_uy_ = static_cast<int>(density_center_y + density_height / 2);
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
  this->setLibName("HYBRID_BOND");
}
Pin *Instance::getHybridBondPin() const {
  return hybrid_bond_pin_;
}
void Instance::setHybridBondPin(Pin *hybrid_bond_pin) {
  if (!hybrid_bond_pin->isHybridBondPin())
    assert(0); // This pin is not a pin for hybrid bond.
  hybrid_bond_pin_ = hybrid_bond_pin;
  connected_pins_.push_back(hybrid_bond_pin);
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
  dbBlock *db_block = db_database_->getChip()->getBlock();
  db_inst_ = dbInst::create(db_block, master, name_.c_str());
}
uint Instance::getCenterX() {
  return getCoordinate().first + floor(getWidth() / 2);
}
uint Instance::getCenterY() {
  return getCoordinate().second + floor(getHeight() / 2);
}
void Instance::assignDie(int die_id) {
  die_id_ = die_id;
}
bool Instance::isInstance() {
  return db_inst_ != nullptr;
}
bool Instance::isDummy() {
  return db_inst_ == nullptr;
}
const vector<Net *> &Instance::getConnectedNets() const {
  return connected_nets_;
}

} // VLSI_backend
