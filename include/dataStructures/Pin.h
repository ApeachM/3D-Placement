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

#ifndef PLACER_INCLUDE_DATASTRUCTURES_PIN_H_
#define PLACER_INCLUDE_DATASTRUCTURES_PIN_H_
#include "db.h"
#include "storages.h"

namespace VLSI_backend {
using namespace std;
class Pin {
 private:
  odb::dbDatabase *parent_ = nullptr;
  odb::dbITerm *db_i_term_ = nullptr;
  odb::dbBTerm *db_b_term_ = nullptr;

  Instance* connected_instance = nullptr;
  Net* connected_net = nullptr;

  bool min_pin_x_field_ = true;
  bool min_pin_y_field_ = true;
  bool max_pin_x_field_ = true;
  bool max_pin_y_field_ = true;

  // density coordinate -- need to refactor
  int offsetCx_ = 0;
  int offsetCy_ = 0;
  int cx_ = 0;
  int cy_ = 0;

  // Please check the equation (4) in the ePlace-MS paper.
  //
  // maxExpSum_: holds exp(x_i/gamma)
  // minExpSum_: holds exp(-x_i/gamma)
  // the x_i is equal to cx_ variable.
  //
  float maxExpSumX_{};
  float maxExpSumY_{};

  float minExpSumX_{};
  float minExpSumY_{};

  // flag variables
  //
  // check whether
  // this pin is considered in a WA models.
  unsigned char hasMaxExpSumX_: 1;
  unsigned char hasMaxExpSumY_: 1;

  unsigned char hasMinExpSumX_: 1;
  unsigned char hasMinExpSumY_: 1;

  bool is_hybrid_bond_pin_ = false;
  Instance* hybrid_bond_ = nullptr;
  Net* intersected_net_ = nullptr;

 public:
  /// Constructors
  Pin() = default;
  explicit Pin(odb::dbITerm *db_iterm);
  explicit Pin(odb::dbBTerm *db_b_term);
  dbITerm *getDbITerm() const;
  dbBTerm *getDbBTerm() const;
  /// return boolean whether it is instance pin or not
  bool isInstancePin();
  /// return boolean whether it is block pin (fixed pad in Die) or not
  bool isBlockPin();

  /// return Instance pointer correspond to the pin
  /// \details
  /// if the pin is block pin, then it returns nullptr.
  Instance *getInstance();

  /// return Net pointer correspond to the pin.
  /// \details
  /// if the signal is POWER or GROUND, then there will be no net (return nullptr).
  Net *getNet();

  /// return the signal type of the pin
  string getSignalType();

  /// return the pin name in the instance
  string getPinName();

  /// return the coordinate of the pin.
  /// \details
  /// the returned coordinate will be the center of the box (pin shape)
  pair<int, int> getCoordinate();;

  Instance *getConnectedInstance() const;
  void setConnectedInstance(Instance *connected_instance);
  Net *getConnectedNet() const;
  void setConnectedNet(Net *connected_net);

  void setMinPinXField(bool min_pin_x_field);
  void setMinPinYField(bool min_pin_y_field);
  void setMaxPinXField(bool max_pin_x_field);
  void setMaxPinYField(bool max_pin_y_field);
  bool isMinPinX() const;
  bool isMinPinY() const;
  bool isMaxPinX() const;
  bool isMaxPinY() const;

  // this function will return the density coordinate
  int cx() { return cx_; }
  int cy() { return cy_; }

  void initDensityCoordinate();
  void updateDensityLocation(Instance *instance);

  // clear WA(Weighted Average) variables.
  void clearWaVars();

  void setMaxExpSumX(float maxExpSumX);
  void setMaxExpSumY(float maxExpSumY);
  void setMinExpSumX(float minExpSumX);
  void setMinExpSumY(float minExpSumY);
  float maxExpSumX() const { return maxExpSumX_; }
  float maxExpSumY() const { return maxExpSumY_; }
  float minExpSumX() const { return minExpSumX_; }
  float minExpSumY() const { return minExpSumY_; }

  bool hasMaxExpSumX() const { return (hasMaxExpSumX_ == 1); }
  bool hasMaxExpSumY() const { return (hasMaxExpSumY_ == 1); }
  bool hasMinExpSumX() const { return (hasMinExpSumX_ == 1); }
  bool hasMinExpSumY() const { return (hasMinExpSumY_ == 1); }

  bool isHybridBondPin() const;
  void setAsHybridBondPin();
  Net *getIntersectedNet() const;
  void setIntersectedNet(Net *intersected_net);
  Instance *getHybridBond() const;
  void setHybridBond(Instance *hybrid_bond);
  void setHybridBondCoordinate(int x, int y);

};
}
#endif //PLACER_INCLUDE_DATASTRUCTURES_PIN_H_
