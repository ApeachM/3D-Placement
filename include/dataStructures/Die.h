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
class Die {
 private:
  dbBlock *db_block_ = nullptr;
  Rect die_shape_{};

  float density_ = 1.0;

  uint width_ = 0;
  uint height_ = 0 ;

  int die_id_ = 0;

 public:
  Die() = default;
  explicit Die(dbBlock *db_block);
  void setDbBlock(dbBlock *db_block);
  void setDieSize(uint width , uint height);
  uint getWidth();
  uint getHeight();

  uint getArea();
  float getDensity() const;
  void setDensity(double density);

  int getLowerLeftX(){
    return die_shape_.ll().getX();
  }
  int getLowerLeftY(){
    return die_shape_.ll().getY();
  }
  int getUpperRightX(){
    return die_shape_.ur().getX();
  }
  int getUpperRightY(){
    return die_shape_.ur().getY();
  }

  int getDieId() const;
  void setDieId(int die_id);

};

} // VLSI_backend

#endif //PLACER_INCLUDE_DATASTRUCTURES_DIE_H_
