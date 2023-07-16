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
#include "storages.h"

namespace flow3D {
using namespace std;
class Instance {
 public:
  /// Constructors
  Instance() = default;
  explicit Instance(odb::dbInst *db_inst);
  Instance(odb::dbInst *db_inst, int id);

  /// return the instance name of the cell
  /// example: _321_
  string getName();

  /// return the library name of the cell
  /// example: DFF_X1
  string getLibName();

  /*!
   * \pre
   * `setLibname()` && `setDbDatabase()`
   * */
  void setLibrary(dbMaster *master);

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
   * this function should be called in only Chip::dataBaseInit() function
   * */
  void setId(int id);

  /// get width of the instance(cell)
  uint getWidth();

  /// get height of the instance(cell)
  uint getHeight();

  /// get area of the instance(cell)
  uint getArea();

  void setWidth(uint width) {
    // this function will be called when making a filler
    width_ = width;
  }
  void setHeight(uint height) {
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
  uint getCenterX();
  uint getCenterY();
  void assignDie(int die_id);

  /*!
   * \brief
   * Determine whether it is dummy cell or not
   * \details
   * If it is not dummy cell, then this returns true.
   * */
  bool isInstance();

  /*!
   * \brief
   * Determine whether it is dummy cell or not
   * \details
   * If it is not dummy cell, then this returns false.
   * */
  bool isDummy();

  // ref: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L139
  bool isFixed();

  void setDensityCenterLocation(int d_cx, int d_cy);

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
  int getDensityCenterX() const { return (d_ux_ + d_lx_) / 2; }
  int getDensityCenterY() const { return (d_uy_ + d_ly_) / 2; }
  int getDensityDeltaX() const { return d_ux_ - d_lx_; }
  int getDensityDeltaY() const { return d_uy_ - d_ly_; }
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
  const vector<Net *> &getConnectedNets() const;
  void setConnectedPins(vector<Pin *> connected_pins);
  void setConnectedNets(vector<Net *> connected_nets);
  void setInstName(const string &name);
  void setLibName(const string &lib_name);
  dbDatabase *getDbDatabase() const;
  void setDbDatabase(dbDatabase *db_database);

 private:
  odb::dbDatabase *db_database_ = nullptr;
  odb::dbInst *db_inst_ = nullptr;

  string name_;
  string libName_;
  int id_{};
  int die_id_ = 0;
  bool is_macro_ = false;
  bool is_locked_ = false;

  vector<Pin *> connected_pins_;
  vector<Net *> connected_nets_;

  /// This is lower left position of instance
  /// This is same with the origin of db_inst_ pointer
  pair<int, int> position_ = pair<int, int>{0, 0};
  uint width_ = 0;
  uint height_ = 0;

  /// density location
  int d_lx_{0};
  int d_ly_{0};
  int d_ux_{0};
  int d_uy_{0};

  // density variables
  float density_scale_{0};
  float gradient_x_{};
  float gradient_y_{};
};
class HybridBond {
 public:
  HybridBond() = default;
  HybridBond(int width, int height, int spacing) {
    width_ = width;
    height_ = height;
    spacing_ = spacing;
  }
  const string &getName() const {
    return name_;
  }
  void setName(const string &name) {
    name_ = name;
  }
  const pair<int, int> &getCoordinate() const {
    return position_;
  }
  void setCoordinate(const pair<int, int> center_position) {
    position_.first = center_position.first - width_ / 2;
    position_.second = center_position.second - height_ / 2;
  }
  Net *getConnectedNet() const {
    return connected_net_;
  }
  void setConnectedNet(Net *connected_net) {
    connected_net_ = connected_net;
  }
  void updatePosition();




 private:
  pair<int, int> position_ = pair<int, int>{0, 0}; // lower left
  Net *connected_net_ = nullptr;
  vector<Pin *> connected_pins_;
  string name_;
  int width_ = 0;
  int height_ = 0;
  int spacing_ = 0;
};
}

#endif //PLACER_INCLUDE_DATASTRUCTURES_INSTANCE_H_























