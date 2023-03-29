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

#ifndef PLACER_INCLUDE_DATASTRUCTURES_DIE_H_
#define PLACER_INCLUDE_DATASTRUCTURES_DIE_H_
#include "db.h"

namespace VLSI_backend {
using namespace odb;
struct RowInfo {
 public:
  int start_x = 0;
  int start_y = 0;
  int row_width = 0;
  int row_height = 0;
  int repeat_count = 0;
};

class Die {
 public:
  Die();
  explicit Die(dbBlock *db_block);
  void setDbBlock(dbBlock *db_block);
  uint getWidth();
  uint getHeight();

  uint getArea();
  float getDensity() const;
  void setDensity(double density);

  int getLowerLeftX() {
    return die_shape_.ll().getX();
  }
  int getLowerLeftY() {
    return die_shape_.ll().getY();
  }
  int getUpperRightX() {
    return die_shape_.ur().getX();
  }
  int getUpperRightY() {
    return die_shape_.ur().getY();
  }

  /*!
   * Caution!! Use this function only when you parsed data without the openDB parse functions.
   * */
  void setDBBasic(string lib_name) {
    assert(db_database_);
    db_database_ = odb::dbDatabase::create();
    db_tech_ = dbTech::create(db_database_);
    db_tech_layer_ = dbTechLayer::create(db_tech_, "Layer", dbTechLayerType::MASTERSLICE);
    db_lib_ = dbLib::create(db_database_, lib_name.c_str(), ',');
    db_chip_ = dbChip::create(db_database_);
    db_block_ = dbBlock::create(db_chip_, (std::to_string(die_id_) + "th Die Block").c_str());

    die_shape_ = db_block_->getDieArea();
    setDieSize(die_shape_.dx(), die_shape_.dy());
  }

  int getDieId() const;
  void setDieId(int die_id);
  dbDatabase *getDbDatabase() const;
  void setDbDatabase(dbDatabase *db_database);
  dbBlock *getDbBlock() const;
  dbTech *getDbTech() const;
  void setDbTech(dbTech *db_tech);
  dbTechLayer *getDbTechLayer() const;
  void setDbTechLayer(dbTechLayer *db_tech_layer);
  dbLib *getDbLib() const;
  void setDbLib(dbLib *db_lib);
  dbChip *getDbChip() const;
  void setDbChip(dbChip *db_chip);
  const string &getTechName() const;
  void setTechName(const string &tech_name);
  int getLibNum() const;
  void setLibNum(int lib_num);
  int getMaxUtil() const;
  void setMaxUtil(int max_util);
  void setRowInfo(int start_x, int start_y, int row_width, int row_height, int repeat_count);

 private:
  dbDatabase *db_database_ = nullptr;
  dbBlock *db_block_ = nullptr;
  dbTech *db_tech_ = nullptr;
  dbTechLayer *db_tech_layer_ = nullptr;
  dbLib *db_lib_ = nullptr;
  dbChip *db_chip_ = nullptr;

  Rect die_shape_{};

  float density_ = 1.0;

  uint width_ = 0;
  uint height_ = 0;

  int die_id_ = 0;
  string tech_name_;
  int lib_num_ = 0;
  int max_util_ = 0; // unit: percent

  RowInfo row_info_;

  void setDieSize(uint width, uint height);
};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_DIE_H_
