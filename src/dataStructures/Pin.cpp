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

#include "Pin.h"

namespace Circuit {

Pin::Pin(odb::dbITerm *db_iterm) {
  db_i_term_ = db_iterm;
  parent_ = db_i_term_->getDb();
}
Pin::Pin(odb::dbBTerm *db_b_term) {
  db_b_term_ = db_b_term;
  parent_ = db_b_term_->getDb();
}
void Pin::setDataStorage(data_storage *data_storage) {
  if (data_storage_ == nullptr)
    data_storage_ = data_storage;
}
void Pin::setDataMapping(data_mapping *data_mapping) {
  if (data_mapping_ == nullptr)
    data_mapping_ = data_mapping;
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
  if (isInstancePin()) {
    return data_mapping_->net_map[db_i_term_->getNet()];
  } else {
    return data_mapping_->net_map[db_b_term_->getNet()];
  }
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
    return data_mapping_->inst_map[db_i_term_->getInst()];
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
} // Circuit