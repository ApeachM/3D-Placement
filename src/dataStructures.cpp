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
#include "NesterovPlacer.h"
#include "Net.h"

namespace flow3D {
// Die //
Die::Die() {
}
Die::Die(dbBlock *db_block) {
  db_block_ = db_block;
  die_shape_ = db_block->getDieArea();
  setDieSize(die_shape_.dx(), die_shape_.dy());
}
void Die::setDbBlock(dbBlock *db_block) {
  db_block_ = db_block;
  die_shape_ = db_block->getDieArea();
  setDieSize(die_shape_.dx(), die_shape_.dy());
}
uint Die::getWidth() {
  return die_shape_.dx();
}
uint Die::getHeight() {
  return die_shape_.dy();
}
uint64 Die::getArea() {
  return die_shape_.area();
}
void Die::setDieSize(uint width, uint height) {
  width_ = width;
  height_ = height;
}
float Die::getDensity() const {
  return density_;
}
void Die::setDensity(double density) {
  density_ = density;
}
int Die::getDieId() const {
  return die_id_;
}
void Die::setDieId(int die_id) {
  die_id_ = die_id;
}
dbBlock *Die::getDbBlock() const {
  return db_block_;
}
dbTech *Die::getDbTech() const {
  return db_tech_;
}
void Die::setDbTech(dbTech *db_tech) {
  db_tech_ = db_tech;
}
dbTechLayer *Die::getDbTechLayer() const {
  return db_tech_layer_;
}
void Die::setDbTechLayer(dbTechLayer *db_tech_layer) {
  db_tech_layer_ = db_tech_layer;
}
dbLib *Die::getDbLib() const {
  return db_lib_;
}
void Die::setDbLib(dbLib *db_lib) {
  db_lib_ = db_lib;
}
dbChip *Die::getDbChip() const {
  return db_chip_;
}
void Die::setDbChip(dbChip *db_chip) {
  db_chip_ = db_chip;
}
const string &Die::getTechName() const {
  return tech_name_;
}
void Die::setTechName(const string &tech_name) {
  tech_name_ = tech_name;
}
dbDatabase *Die::getDbDatabase() const {
  return db_database_;
}
void Die::setDbDatabase(dbDatabase *db_database) {
  db_database_ = db_database;
}
int Die::getLibNum() const {
  return lib_num_;
}
void Die::setLibNum(int lib_num) {
  lib_num_ = lib_num;
}
int Die::getMaxUtil() const {
  return max_util_;
}
void Die::setMaxUtil(int max_util) {
  max_util_ = max_util;
}
void Die::setRowInfo(int start_x, int start_y, int row_width, int row_height, int repeat_count) {
  row_info_.start_x = start_x;
  row_info_.start_y = start_y;
  row_info_.row_width = row_width;
  row_info_.row_height = row_height;
  row_info_.repeat_count = repeat_count;
}
const RowInfo &Die::getRowInfo() const {
  return row_info_;
}

// Instance //
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
void HybridBond::updatePosition() {
  int center_x = static_cast<int>((connected_net_->ux() + connected_net_->lx()) / 2);
  int center_y = static_cast<int>((connected_net_->uy() + connected_net_->ly()) / 2);
  setCoordinate({center_x, center_y});
}

// Net //
Net::Net(odb::dbNet *db_net) {
  db_database_ = db_net->getDb();
  db_net_ = db_net;
  name_ = db_net->getName();
}
dbNet *Net::getDbNet() const {
  if (db_net_ == nullptr)
    std::cout << "Invalid access to dbNet from net" << std::endl;
  return db_net_;
}
vector<Pin *> Net::getConnectedPins() {
  // TODO: can be more simplified?
  vector<Pin *> pins;
  for (Pin *pin : connected_pins_) {
    pins.push_back(pin);
  }
  return pins;
}
int Net::getWeight() {
  return db_net_->getWeight();
}
string Net::getSignalType() {
  return db_net_->getSigType().getString();
}
string Net::getName() {
  return db_net_->getName();
}
ulong Net::getHPWL() {
  Rect net_box;
  net_box.mergeInit();

  for (dbITerm *iterm : db_net_->getITerms()) {
    int x, y;
    if (iterm->getAvgXY(&x, &y)) {
      Rect iterm_rect(x, y, x, y);
      net_box.merge(iterm_rect);
    } else {
      // This clause is sort of worthless because getAvgXY prints
      // a warning when it fails.
      dbInst *inst = iterm->getInst();
      dbBox *inst_box = inst->getBBox();
      int center_x = (inst_box->xMin() + inst_box->xMax()) / 2;
      int center_y = (inst_box->yMin() + inst_box->yMax()) / 2;
      Rect inst_center(center_x, center_y, center_x, center_y);
      net_box.merge(inst_center);
    }
  }

  for (dbBTerm *bterm : db_net_->getBTerms()) {
    for (dbBPin *bpin : bterm->getBPins()) {
      // TODO: Debug here
      dbPlacementStatus status = bpin->getPlacementStatus();
      if (status.isPlaced()) {
        Rect pin_bbox = bpin->getBBox();
        int center_x = (pin_bbox.xMin() + pin_bbox.xMax()) / 2;
        int center_y = (pin_bbox.yMin() + pin_bbox.yMax()) / 2;
        Rect pin_center(center_x, center_y, center_x, center_y);
        net_box.merge(pin_center);
      }
    }
  }

  return net_box.dx() + net_box.dy();
}
ulong Net::getHPWL(DIE_ID die_id) {
  Rect net_box;
  net_box.mergeInit();

  for (Pin *pin : this->getConnectedPins()) {
    int x, y;
    if (pin->isInstancePin()) {
      if (pin->getInstance()->getDieId() != die_id)
        continue;
      dbITerm *i_term = pin->getDbITerm();
      assert(i_term);
      if (i_term->getAvgXY(&x, &y)) {
        Rect i_term_rect(x, y, x, y);
        net_box.merge(i_term_rect);
      } else {
        // This clause is sort of worthless because getAvgXY prints
        // a warning when it fails.
        dbInst *inst = i_term->getInst();
        dbBox *inst_box = inst->getBBox();
        int center_x = (inst_box->xMin() + inst_box->xMax()) / 2;
        int center_y = (inst_box->yMin() + inst_box->yMax()) / 2;
        Rect inst_center(center_x, center_y, center_x, center_y);
        net_box.merge(inst_center);
      }
    }
//    else if (pin->isBlockPin()) {
//      dbBTerm *b_term = pin->getDbBTerm();
//      assert(b_term);
//      for (dbBPin *b_pin : b_term->getBPins()) {
//        dbPlacementStatus status = b_pin->getPlacementStatus();
//        if (status.isPlaced()) {
//          Rect pin_bbox = b_pin->getBBox();
//          int center_x = (pin_bbox.xMin() + pin_bbox.xMax()) / 2;
//          int center_y = (pin_bbox.yMin() + pin_bbox.yMax()) / 2;
//          Rect pin_center(center_x, center_y, center_x, center_y);
//          net_box.merge(pin_center);
//        }
//      }
//    } else {
//      assert(0);
//    }
  }

  pair<int, int> hb_coordinate = hybrid_bond_->getCoordinate();
  Rect hybrid_bond_rect(hb_coordinate.first, hb_coordinate.second, hb_coordinate.first, hb_coordinate.second);
  net_box.merge(hybrid_bond_rect);

  return net_box.dx() + net_box.dy();
}
int64_t Net::hpwl() {
  return static_cast<int64_t>((ux_ - lx_) + (uy_ - ly_));
}
// TODO: check is that right or not
void Net::updateBox(int die_ID, bool consider_other_die) {
  lx_ = ly_ = INT_MAX;
  ux_ = uy_ = INT_MIN;
  if (consider_other_die || die_ID == 0) {
    for (Pin *pin : getConnectedPins()) {
      lx_ = std::min(pin->cx(), lx_);
      ly_ = std::min(pin->cy(), ly_);
      ux_ = std::max(pin->cx(), ux_);
      uy_ = std::max(pin->cy(), uy_);
    }
  } else {
    // consider only for the one die
    // ignore the pin for this net for updating box
    for (Pin *pin : getConnectedPins()) {
      if (pin->isBlockPin()) {
        lx_ = std::min(pin->cx(), lx_);
        ly_ = std::min(pin->cy(), ly_);
        ux_ = std::max(pin->cx(), ux_);
        uy_ = std::max(pin->cy(), uy_);
      } else if (die_ID == pin->getInstance()->getDieId()) {
        // can not access the instance id when the pin is block pin.
        lx_ = std::min(pin->cx(), lx_);
        ly_ = std::min(pin->cy(), ly_);
        ux_ = std::max(pin->cx(), ux_);
        uy_ = std::max(pin->cy(), uy_);
      }
    }
  }
  assert(lx_ <= ux_);
  assert(ly_ <= uy_);
}
bool Net::isIntersected() const {
  return intersected_;
}
void Net::setAsIntersected() {
  Net::intersected_ = true;
  die_id_ = DIE_ID::INTERSECTED;
}
int Net::getDieId() const {
  return die_id_;
}
void Net::setDieId(int die_id) {
  die_id_ = die_id;
}
HybridBond *Net::getHybridBond() const {
  if (intersected_ == false)
    assert(0);
  return hybrid_bond_;
}
void Net::setHybridBond(HybridBond *hybrid_bond_pin) {
  assert(intersected_ == true);
  hybrid_bond_ = hybrid_bond_pin;
}
void Net::setConnectedPins(const vector<Pin *> &connected_pins) {
  connected_pins_ = connected_pins;
}
const vector<Instance *> &Net::getConnectedInstances() const {
  return connected_instances_;
}
void Net::setConnectedInstances(const vector<Instance *> &connected_instances) {
  connected_instances_ = connected_instances;
}
pair<int, int> Net::getCenter() {
  // This function is similar to the updateBox, but that function requires calling many functions of the Nesterov placer
  // For being independent of that, I made this function temporary.
  // The net variables for Nesterov should be separated considering the algorithms.
  // We should do that, but pending it in the TODO.

  int lower_left_x, lower_left_y;
  int upper_right_x, upper_right_y;
  lower_left_x = lower_left_y = INT_MAX;
  upper_right_x = upper_right_y = INT_MIN;

  for (Pin *pin : getConnectedPins()) {
    pair<int, int> coordinate = pin->getCoordinate();
    lower_left_x = std::min(coordinate.first, lower_left_x);
    lower_left_y = std::min(coordinate.second, lower_left_y);
    upper_right_x = std::max(coordinate.first, upper_right_x);
    upper_right_y = std::max(coordinate.second, upper_right_y);
  }
  assert(lower_left_x <= upper_right_x);
  assert(lower_left_y <= upper_right_y);

  int center_x = floor((lower_left_x + upper_right_x) / 2);
  int center_y = floor((lower_left_y + upper_right_y) / 2);

  return {center_x, center_y};
}

// Pin //
Pin::Pin(odb::dbITerm *db_iterm) {
  db_i_term_ = db_iterm;
  parent_ = db_i_term_->getDb();
}
Pin::Pin(odb::dbBTerm *db_b_term) {
  db_b_term_ = db_b_term;
  parent_ = db_b_term_->getDb();
}
dbITerm *Pin::getDbITerm() const {
  return db_i_term_;
}
dbBTerm *Pin::getDbBTerm() const {
  return db_b_term_;
}
bool Pin::isInstancePin() {
  if (db_i_term_ != nullptr)
    return true;
  else
    return false;
}
bool Pin::isBlockPin() {
  if (db_b_term_ != nullptr)
    return true;
  else
    return false;
}
Net *Pin::getNet() {
  return connected_net;
}
string Pin::getSignalType() {
  if (isInstancePin())
    return db_i_term_->getSigType().getString();
  else
    return db_b_term_->getSigType().getString();
}
pair<int, int> Pin::getCoordinate() {
  int x = 0, y = 0;
  if (isInstancePin()) {
    db_i_term_->getAvgXY(&x, &y);
  } else if (isBlockPin()) {
    if (db_b_term_->getBPins().size() > 1)
      assert(1);  // if this case occurs, please contact to me :(
    for (auto bPin : db_b_term_->getBPins()) {
      if (bPin->getBoxes().size() > 1)
        assert(1);  // if this case occurs, please contact to me :(
      for (auto box : bPin->getBoxes()) {
        x = (box->xMin() + (int) box->getDX());
        y = (box->yMin() + (int) box->getDY());
      }
    }
  }
  return pair<int, int>{x, y};
}
Instance *Pin::getInstance() {
  if (isInstancePin())
    return connected_instance;
  else
    return nullptr;
}
string Pin::getPinName() {
  string name;
  if (isInstancePin()) {
    name = getDbITerm()->getMTerm()->getName();
  } else if (isBlockPin()) {
    name = getDbBTerm()->getName();
  } else {
    assert(0);
  }
  return name;
}
void Pin::initDensityCoordinate() {
  cx_ = getCoordinate().first;
  cy_ = getCoordinate().second;
  if (db_i_term_) {
    offsetCx_ = cx_ - getInstance()->getCoordinate().first;
    offsetCy_ = cy_ - getInstance()->getCoordinate().second;
  } else if (db_b_term_) {
    // TODO: ??? is this right ???
    offsetCx_ = cx_;
    offsetCy_ = cy_;
  } else {
    // if this is a pin for hybrid
    offsetCx_ = cx_;
    offsetCy_ = cy_;
  }
}
void Pin::setMinPinXField(bool min_pin_x_field) {
  min_pin_x_field_ = min_pin_x_field;
}
void Pin::setMinPinYField(bool min_pin_y_field) {
  min_pin_y_field_ = min_pin_y_field;
}
void Pin::setMaxPinXField(bool max_pin_x_field) {
  max_pin_x_field_ = max_pin_x_field;
}
void Pin::setMaxPinYField(bool max_pin_y_field) {
  max_pin_y_field_ = max_pin_y_field;
}
bool Pin::isMinPinX() const {
  return min_pin_x_field_;
}
bool Pin::isMinPinY() const {
  return min_pin_y_field_;
}
bool Pin::isMaxPinX() const {
  return max_pin_x_field_;
}
bool Pin::isMaxPinY() const {
  return max_pin_y_field_;
}
Instance *Pin::getConnectedInstance() const {
  return connected_instance;
}
void Pin::setConnectedInstance(Instance *connected_instance) {
  Pin::connected_instance = connected_instance;
}
Net *Pin::getConnectedNet() const {
  return connected_net;
}
void Pin::setConnectedNet(Net *connected_net) {
  Pin::connected_net = connected_net;
}
void Pin::updateDensityLocation(Instance *instance) {
  // why is this "instance->getDensityCenterX + offsetCx"?
  // shouldn't it be "instance->dLy + offsetCx"?
  cx_ = instance->getDensityCenterX() + offsetCx_;
  cy_ = instance->getDensityCenterY() + offsetCy_;
}
void Pin::clearWaVars() {
  hasMaxExpSumX_ = 0;
  hasMaxExpSumY_ = 0;
  hasMinExpSumX_ = 0;
  hasMinExpSumY_ = 0;

  maxExpSumX_ = maxExpSumY_ = 0;
  minExpSumX_ = minExpSumY_ = 0;
}
void Pin::setMaxExpSumX(float maxExpSumX) {
  hasMaxExpSumX_ = 1;
  maxExpSumX_ = maxExpSumX;
  assert(!(isnan(maxExpSumX) || isinf(maxExpSumX)));
}
void Pin::setMaxExpSumY(float maxExpSumY) {
  hasMaxExpSumY_ = 1;
  maxExpSumY_ = maxExpSumY;
  assert(!(isnan(maxExpSumY) || isinf(maxExpSumY)));
}
void Pin::setMinExpSumX(float minExpSumX) {
  hasMinExpSumX_ = 1;
  minExpSumX_ = minExpSumX;
  assert(!(isnan(minExpSumX) || isinf(minExpSumX)));
}
void Pin::setMinExpSumY(float minExpSumY) {
  hasMinExpSumY_ = 1;
  minExpSumY_ = minExpSumY;
  assert(!(isnan(minExpSumY) || isinf(minExpSumY)));
}

// Bin //
Chip::NesterovPlacer::Bin::Bin()
    : x_(0),
      y_(0),
      lx_(0),
      ly_(0),
      ux_(0),
      uy_(0),
      nonPlaceArea_(0),
      instPlacedArea_(0),
      instPlacedAreaUnscaled_(0),
      nonPlaceAreaUnscaled_(0),
      fillerArea_(0),
      density_(0),
      targetDensity_(0),
      electroPhi_(0),
      electroForceX_(0),
      electroForceY_(0) {
}
Chip::NesterovPlacer::Bin::Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity) : Bin() {
  x_ = x;
  y_ = y;
  lx_ = lx;
  ly_ = ly;
  ux_ = ux;
  uy_ = uy;
  targetDensity_ = targetDensity;
}
Chip::NesterovPlacer::Bin::~Bin() {
  x_ = y_ = 0;
  lx_ = ly_ = ux_ = uy_ = 0;
  nonPlaceArea_ = instPlacedArea_ = fillerArea_ = nonPlaceAreaUnscaled_
      = instPlacedAreaUnscaled_ = 0;
  electroPhi_ = electroForceX_ = electroForceY_ = 0;
  density_ = targetDensity_ = 0;
}
const int64_t Chip::NesterovPlacer::Bin::binArea() const {
  return static_cast<int64_t>(dx()) * static_cast<int64_t>(dy());
}
float Chip::NesterovPlacer::Bin::density() const {
  return density_;
}
float Chip::NesterovPlacer::Bin::targetDensity() const {
  return targetDensity_;
}
float Chip::NesterovPlacer::Bin::electroForceX() const {
  return electroForceX_;
}
float Chip::NesterovPlacer::Bin::electroForceY() const {
  return electroForceY_;
}
float Chip::NesterovPlacer::Bin::electroPhi() const {
  return electroPhi_;
}
void Chip::NesterovPlacer::Bin::setDensity(float density) {
  density_ = density;
}
void Chip::NesterovPlacer::Bin::setTargetDensity(float density) {
  targetDensity_ = density;
}
void Chip::NesterovPlacer::Bin::setElectroForce(float electroForceX, float electroForceY) {
  electroForceX_ = electroForceX;
  electroForceY_ = electroForceY;
}
void Chip::NesterovPlacer::Bin::setElectroPhi(float phi) {
  electroPhi_ = phi;
}
int Chip::NesterovPlacer::Bin::x() const {
  return x_;
}
int Chip::NesterovPlacer::Bin::y() const {
  return y_;
}
int Chip::NesterovPlacer::Bin::lx() const {
  return lx_;
}
int Chip::NesterovPlacer::Bin::ly() const {
  return ly_;
}
int Chip::NesterovPlacer::Bin::ux() const {
  return ux_;
}
int Chip::NesterovPlacer::Bin::uy() const {
  return uy_;
}
int Chip::NesterovPlacer::Bin::cx() const {
  return (ux_ + lx_) / 2;
}
int Chip::NesterovPlacer::Bin::cy() const {
  return (uy_ + ly_) / 2;
}
int Chip::NesterovPlacer::Bin::dx() const {
  return (ux_ - lx_);
}
int Chip::NesterovPlacer::Bin::dy() const {
  return (uy_ - ly_);
}
void Chip::NesterovPlacer::Bin::setNonPlaceArea(int64_t area) {
  nonPlaceArea_ = area;
}
void Chip::NesterovPlacer::Bin::setNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ = area;
}
void Chip::NesterovPlacer::Bin::setInstPlacedArea(int64_t area) {
  instPlacedArea_ = area;
}
void Chip::NesterovPlacer::Bin::setInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ = area;
}
void Chip::NesterovPlacer::Bin::setFillerArea(int64_t area) {
  fillerArea_ = area;
}
void Chip::NesterovPlacer::Bin::addNonPlaceArea(int64_t area) {
  nonPlaceArea_ += area;
}
void Chip::NesterovPlacer::Bin::addInstPlacedArea(int64_t area) {
  instPlacedArea_ += area;
}
void Chip::NesterovPlacer::Bin::addNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ += area;
}
void Chip::NesterovPlacer::Bin::addInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ += area;
}
void Chip::NesterovPlacer::Bin::addFillerArea(int64_t area) {
  fillerArea_ += area;
}
}