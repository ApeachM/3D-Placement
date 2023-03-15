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
#include "Die.h"

namespace VLSI_backend {
Die::Die() {
}
Die::Die(dbBlock *db_block) {
  db_block_ = db_block;
  db_block->getDieArea(die_shape_);
  setDieSize(die_shape_.dx(), die_shape_.dy());
}
void Die::setDbBlock(dbBlock *db_block) {
  db_block_ = db_block;
  db_block->getDieArea(die_shape_);
  setDieSize(die_shape_.dx(), die_shape_.dy());
}
uint Die::getWidth() {
  return die_shape_.dx();
}
uint Die::getHeight() {
  return die_shape_.dy();
}
uint Die::getArea() {
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
} // VLSI_backend