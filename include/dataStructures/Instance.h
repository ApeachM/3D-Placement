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

namespace Circuit {
using namespace std;
class Instance {
 private:
  odb::dbDatabase *db_database_ = nullptr;
  odb::dbInst *db_inst_ = nullptr;
  data_storage *data_storage_ = nullptr;
  data_mapping *data_mapping_ = nullptr;

  string name_;
  string libName_;

  /// This is lower left position of instance
  /// This is same with the origin of db_inst_ pointer
  pair<int, int> position_ = pair<int, int>{0, 0};
 public:

  /// Constructors
  Instance() = default;
  explicit Instance(odb::dbInst *db_inst);
  Instance(odb::dbInst *db_inst, data_storage *data_storage, data_mapping *data_mapping);

  /// set data mapping pointer
  void setDataMapping(data_mapping *data_mapping);

  /// set data storage pointer
  void setDataStorage(data_storage *data_storage);

  /// return the instance name of the cell
  /// example: _321_
  string getName();

  /// return the library name of the cell
  /// example: DFF_X1
  string getLibName();

  /// get width of the instance(cell)
  uint getWidth();

  /// get height of the instance(cell)
  uint getHeight();

  /// get area of the instance(cell)
  uint getArea();

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


};

}

#endif //PLACER_INCLUDE_DATASTRUCTURES_INSTANCE_H_























