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

namespace VLSI_backend {

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
void Net::setDataMapping(data_mapping *data_mapping) {
  if (data_mapping_ == nullptr)
    data_mapping_ = data_mapping;
}
void Net::setDataStorage(data_storage *data_storage) {
  if (data_storage_ == nullptr)
    data_storage_ = data_storage;
}
vector<Pin *> Net::getConnectedPins() {
  vector<Pin *> pins;
  for (dbITerm *db_i_term : db_net_->getITerms()) {
    pins.push_back(data_mapping_->pin_map_i[db_i_term]);
  }
  for (dbBTerm *db_b_term : db_net_->getBTerms()) {
    pins.push_back(data_mapping_->pin_map_b[db_b_term]);
  }
  if (isIntersected()) {
    pins.push_back(hybrid_bond_pin_);
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
int64_t Net::hpwl() {
  return static_cast<int64_t>((ux_ - lx_) + (uy_ - ly_));
}
// TODO: check is that right or not
void Net::updateBox(int die_ID) {
  lx_ = ly_ = INT_MAX;
  ux_ = uy_ = INT_MIN;
  // TODO: check the default value is 0 or not
  if (!isIntersected()) {
    for (Pin *gPin : getConnectedPins()) {
      lx_ = std::min(gPin->cx(), lx_);
      ly_ = std::min(gPin->cy(), ly_);
      ux_ = std::max(gPin->cx(), ux_);
      uy_ = std::max(gPin->cy(), uy_);
    }
  } else {
    for (Pin *gPin : getConnectedPins()) {
      if (!gPin->isBlockPin()) {
        if (die_ID == gPin->getInstance()->getDieId()) {
          lx_ = std::min(gPin->cx(), lx_);
          ly_ = std::min(gPin->cy(), ly_);
          ux_ = std::max(gPin->cx(), ux_);
          uy_ = std::max(gPin->cy(), uy_);
        }
      } else {
        // TODO: need to more accuracy method for block pins
        lx_ = std::min(gPin->cx(), lx_);
        ly_ = std::min(gPin->cy(), ly_);
        ux_ = std::max(gPin->cx(), ux_);
        uy_ = std::max(gPin->cy(), uy_);
      }
    }
  }
}
bool Net::isIntersected() const {
  return intersected_;
}
void Net::setAsIntersected() {
  Net::intersected_ = true;
  die_id_ = -1;
}
int Net::getDieId() const {
  return die_id_;
}
void Net::setDieId(int die_id) {
  die_id_ = die_id;
}
Pin *Net::getHybridBondPin() const {
  if (intersected_ == false)
    assert(0);
  return hybrid_bond_pin_;
}
void Net::setHybridBondPin(Pin *hybrid_bond_pin) {
  if (intersected_ == false)
    assert(0);
  hybrid_bond_pin_ = hybrid_bond_pin;
}
} // VLSI_backend