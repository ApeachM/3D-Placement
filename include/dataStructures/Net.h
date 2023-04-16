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

#ifndef PLACER_INCLUDE_DATASTRUCTURES_NET_H_
#define PLACER_INCLUDE_DATASTRUCTURES_NET_H_

#include <string>
#include <iostream>
#include "db.h"
#include "storages.h"
namespace VLSI_backend {
using namespace std;
class Net {
 public:
  /// Constructors
  Net() = default;
  explicit Net(odb::dbNet *db_net);
  /// methods for Chip.init()
  odb::dbNet *getDbNet() const;

  /// get net name
  string getName();

  /// get the connected pin pointers in the net
  vector<Pin *> getConnectedPins();

  /// get signal type of the net
  /// \example
  /// \c SIGNAL, \c POWER, or \c GROUND
  string getSignalType();

  /// get weight of the net
  int getWeight();

  /// get HPWLe of the net
  ulong getHPWL();

  inline void addWaExpMinSumX(float waExpMinSumX) { waExpMinSumX_ += waExpMinSumX; }
  inline void addWaXExpMinSumX(float waXExpMinSumX) { waXExpMinSumX_ += waXExpMinSumX; }
  inline void addWaExpMinSumY(float waExpMinSumY) { waExpMinSumY_ += waExpMinSumY; }
  inline void addWaYExpMinSumY(float waExpXMinSumY) { waYExpMinSumY_ += waExpXMinSumY; }
  inline void addWaExpMaxSumX(float waExpMaxSumX) { waExpMaxSumX_ += waExpMaxSumX; }
  inline void addWaXExpMaxSumX(float waXExpMaxSumX) { waXExpMaxSumX_ += waXExpMaxSumX; }
  inline void addWaExpMaxSumY(float waExpMaxSumY) { waExpMaxSumY_ += waExpMaxSumY; }
  inline void addWaYExpMaxSumY(float waYExpMaxSumY) { waYExpMaxSumY_ += waYExpMaxSumY; }
  inline float waExpMinSumX() const { return waExpMinSumX_; }
  inline float waXExpMinSumX() const { return waXExpMinSumX_; }
  inline float waExpMinSumY() const { return waExpMinSumY_; }
  inline float waYExpMinSumY() const { return waYExpMinSumY_; }
  inline float waExpMaxSumX() const { return waExpMaxSumX_; }
  inline float waXExpMaxSumX() const { return waXExpMaxSumX_; }
  inline float waExpMaxSumY() const { return waExpMaxSumY_; }
  inline float waYExpMaxSumY() const { return waYExpMaxSumY_; }
  void clearWaVars() {
    waExpMinSumX_ = 0;
    waXExpMinSumX_ = 0;

    waExpMaxSumX_ = 0;
    waXExpMaxSumX_ = 0;

    waExpMinSumY_ = 0;
    waYExpMinSumY_ = 0;

    waExpMaxSumY_ = 0;
    waYExpMaxSumY_ = 0;
  }
  inline int lx() const { return lx_; }
  inline int ly() const { return ly_; }
  inline int ux() const { return ux_; }
  inline int uy() const { return uy_; }
  void updateBox(int die_ID = 0);
  float totalWeight() const { return timingWeight_ * customWeight_; }
  float timingWeight() const { return timingWeight_; }
  float customWeight() const { return customWeight_; }
  int64_t hpwl();

  bool isIntersected() const;
  void setAsIntersected();
  int getDieId() const;
  void setDieId(int die_id);
  Pin *getHybridBondPin() const;
  void setHybridBondPin(Pin *hybrid_bond_pin);

  void setConnectedPins(const vector<Pin *> &connected_pins);
  const vector<Instance *> &getConnectedInstances() const;
  void setConnectedInstances(const vector<Instance *> &connected_instances);

 private:
  odb::dbDatabase *db_database_ = nullptr;
  odb::dbNet *db_net_ = nullptr;

  vector<Pin*> connected_pins_;
  vector<Instance*> connected_instances_;

  std::string name_;

  int lx_;
  int ly_;
  int ux_;
  int uy_;

  float timingWeight_ = 1;
  float customWeight_ = 1;

  /*
   weighted average WL model stor for better indexing
   Please check the equation (4) in the ePlace-MS paper.

   WA: weighted Average
   saving four variable will be helpful for
   calculating the WA gradients/wirelengths.

   gamma: modeling accuracy.
  */

  /*
   X forces.
   waExpMinSumX_: store sigma {exp(x_i/gamma)}
   waXExpMinSumX_: store signa {x_i*exp(e_i/gamma)}
   waExpMaxSumX_ : store sigma {exp(-x_i/gamma)}
   waXExpMaxSumX_: store sigma {x_i*exp(-x_i/gamma)}
  */
  float waExpMinSumX_ = 0;
  float waXExpMinSumX_ = 0;

  float waExpMaxSumX_ = 0;
  float waXExpMaxSumX_ = 0;

  /*
   Y forces.
   waExpMinSumY_: store sigma {exp(y_i/gamma)}
   waYExpMinSumY_: store signa {y_i*exp(e_i/gamma)}
   waExpMaxSumY_ : store sigma {exp(-y_i/gamma)}
   waYExpMaxSumY_: store sigma {y_i*exp(-y_i/gamma)}
  */
  float waExpMinSumY_ = 0;
  float waYExpMinSumY_ = 0;

  float waExpMaxSumY_ = 0;
  float waYExpMaxSumY_ = 0;

  unsigned char isDontCare_: 1;

  // if this net has been shared by two die, then we pretend it is intersected_
  bool intersected_ = false;
  // if this net is on the die1 or die2, then this value will be 1 or 2, respectively.
  // if this net is intersected, then this value will be -1.
  // if this net is on virtual, then this value will be zero.
  int die_id_ = 0;
  // this should be nullptr only when `intersected_` is false
  Pin* hybrid_bond_pin_ = nullptr;

};

}

#endif //PLACER_INCLUDE_DATASTRUCTURES_NET_H_
