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
#include "structures.h"

namespace Circuit {
using namespace std;
class Pin {
 private:
  odb::dbDatabase *parent_ = nullptr;
  odb::dbITerm *db_i_term_ = nullptr;
  odb::dbBTerm *db_b_term_ = nullptr;

  data_storage *data_storage_ = nullptr;
  data_mapping *data_mapping_ = nullptr;

 public:
  /// Constructors
  Pin() = default;
  explicit Pin(odb::dbITerm *db_iterm);
  explicit Pin(odb::dbBTerm *db_b_term);
  /// methods for Circuit.init()
  void setDataStorage(data_storage *data_storage);
  void setDataMapping(data_mapping *data_mapping);
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

};
}
#endif //PLACER_INCLUDE_DATASTRUCTURES_PIN_H_
