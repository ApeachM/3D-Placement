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

Die::Die(dbBlock *db_block) {
  db_block_ = db_block;
  db_block->getDieArea(die_shape_);
  setDieSize( die_shape_.dx(), die_shape_.dy());
}
void Die::setDbBlock(dbBlock *db_block) {
  db_block_ = db_block;
  db_block->getDieArea(die_shape_);
  setDieSize( die_shape_.dx(), die_shape_.dy());
}
uint Die::getWidth() {
  return  die_shape_.dx();
}
uint Die::getHeight() {
  return  die_shape_.dy();
}
uint Die::getArea() {
  return die_shape_.area();
}
void Die::setDieSize(uint width, uint height) {
  // TODO: set the size of dbBlock, and call setDbBlock
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
} // VLSI_backend