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
#include "structures.h"
namespace VLSI_backend {
using namespace std;
class Pin;
class Instance;
class Net {
 private:
  odb::dbDatabase *db_database_ = nullptr;
  odb::dbNet *db_net_ = nullptr;
  data_storage *data_storage_ = nullptr;
  data_mapping *data_mapping_ = nullptr;

  std::string name_;

 public:
  /// Constructors
  Net() = default;
  explicit Net(odb::dbNet *db_net);
  /// methods for Circuit.init()
  dbNet *getDbNet() const;
  void setDataMapping(data_mapping *data_mapping);
  void setDataStorage(data_storage *data_storage);

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

};

}

#endif //PLACER_INCLUDE_DATASTRUCTURES_NET_H_
