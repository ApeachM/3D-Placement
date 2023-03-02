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


#ifndef PLACER_INCLUDE_DATASTRUCTURES_INSTANCE_H_
#define PLACER_INCLUDE_DATASTRUCTURES_INSTANCE_H_
#include <iostream>
#include <string>
#include "db.h"
#include "structures.h"

namespace VLSI_backend {
using namespace std;
class Instance {
 private:
  odb::dbDatabase *db_database_ = nullptr;
  odb::dbInst *db_inst_ = nullptr;

  string name_;
  string libName_;
  int id_{};
  int die_id_ = 0;
  bool is_macro_ = false;
  bool is_locked_ = false;
  bool is_hybrid_bond_ = false;
  Pin* hybrid_bond_pin_ = nullptr;

  vector<Pin*> connected_pins_;
  vector<Net*> connected_nets_;

  /// This is lower left position of instance
  /// This is same with the origin of db_inst_ pointer
  pair<int, int> position_ = pair<int, int>{0, 0};
  uint width_ = 0;
  uint height_ = 0;

  /// density location
  int dLx_{0};
  int dLy_{0};
  int dUx_{0};
  int dUy_{0};

  // density variables
  float densityScale_{0};
  float gradientX_{};
  float gradientY_{};

 public:

  /// Constructors
  Instance() = default;
  explicit Instance(odb::dbInst *db_inst);

  /// return the instance name of the cell
  /// example: _321_
  string getName();

  /// return the library name of the cell
  /// example: DFF_X1
  string getLibName();

  /*!
   * \brief
   * return the cell id.
   * example: 13
   * */
  int getId() const;

  /*!
   * \brief
   * set the cell id
   * \details
   * this function should be called in only Chip::init() function
   * */
  void setId(int id);

  /// get width of the instance(cell)
  uint getWidth();

  /// get height of the instance(cell)
  uint getHeight();

  /// get area of the instance(cell)
  uint getArea();

  void setWidth(uint width){
    // this function will be called when making a filler
    width_ = width;
  }
  void setHeight(uint height){
    // this function will be called when making a filler
    height_ = height;
  }

  std::vector<Pin *> getPins();

  /// get db instance pointer
  dbInst *getDbInst() const;

  /// get coordinate of the instance(cell)
  /// \details
  /// this function will return the lower left coordinate of instance.
  pair<int, int> getCoordinate();

  /// set coordinate of the instance(cell)
  /// \details
  /// this function will set the coordinate of instance as int data type.
  /// after calling this function, the pins correspond to the cell will be also moved automatically.
  void setCoordinate(int x, int y);
  /// check whether it is placed or not
  bool isPlaced();

  uint getCenterX() {
    return getCoordinate().first + floor(getWidth() / 2);
  }
  uint getCenterY() {
    return getCoordinate().second + floor(getHeight() / 2);
  }

  void assignDie(int die_id) {
    die_id_ = die_id;
  }

  /*!
   * \brief
   * Determine whether it is dummy cell or not
   * \details
   * If it is not dummy cell, then this returns true.
   * */
  bool isInstance() {
    return db_inst_ != nullptr;
  }

  /*!
   * \brief
   * Determine whether it is dummy cell or not
   * \details
   * If it is not dummy cell, then this returns false.
   * */
  bool isDummy() {
    return db_inst_ == nullptr;
  }

  // ref: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L139
  bool isFixed();

  void setDensityCenterLocation(int dCx, int dCy);

  /// getter and setters for density variables
  int dLx() const;
  void setDLx(int d_lx);
  int dLy() const;
  void setDLy(int d_ly);
  int dUx() const;
  void setDUx(int d_ux);
  int dUy() const;
  void setDUy(int d_uy);
  float densityScale() const;
  void setDensityScale(float density_scale);
  float getGradientX() const;
  void setGradientX(float gradient_x);
  float getGradientY() const;
  void setGradientY(float gradient_y);
  int getDensityCenterX() const { return (dUx_ + dLx_) / 2; }
  int getDensityCenterY() const { return (dUy_ + dLy_) / 2; }
  int getDensityDeltaX() const { return dUx_ - dLx_; }
  int getDensityDeltaY() const { return dUy_ - dLy_; }
  int lx() { return getCoordinate().first; }
  int ly() { return getCoordinate().second; }
  int ux() { return getCoordinate().first + getWidth(); }
  int uy() { return getCoordinate().second + getHeight(); }
  int cx() { return lx() + getWidth() / 2; }
  int cy() { return ly() + getHeight() / 2; }
  int dx() { return ux() - lx(); }
  int dy() { return uy() - ly(); }
  void setDensityValueAsDefault();
  bool isMacro() const { return is_macro_; }
  bool isLocked() const { return is_locked_; }
  void setDensitySize(float density_width, float density_height);
  bool isMacroInstance();
  bool isStdInstance();
  bool isFiller();
  void setDensityLocation(float dLx, float dLy);
  int getDieId() const;
  bool isHybridBond() const;
  void setAsHybridBond();
  Pin *getHybridBondPin() const;
  void setHybridBondPin(Pin *hybrid_bond_pin);
  void setConnectedPins(vector<Pin *> connected_pins);
  void setConnectedNets(vector<Net *> connected_nets);
};

}

#endif //PLACER_INCLUDE_DATASTRUCTURES_INSTANCE_H_























