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
} // VLSI_backend