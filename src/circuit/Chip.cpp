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

using namespace std;

namespace VLSI_backend {
void Chip::parse(const string &lef_name, const string &def_name) {
  parser_.readLef(lef_name);
  parser_.readDef(def_name);
  this->init();
}
void Chip::init() {
  dbBlock *block = parser_.db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();
  dbSet<dbNet> db_nets = block->getNets();

  /*!
   * @brief
   * Instance setting
   *
   * @details
   * 1. It makes real instance data and store in \c data_storage.instances. \n
   * 2. Then it makes pointer set for \c Chip class, \n
   * 3. And also makes mapping from \c db_instance to instance pointer. \n\n
   * */
  data_storage_.instances.reserve(db_instances.size());  // real data for instance
  instance_pointers_.reserve(db_instances.size());  // pointer data for instances
  // 1. make real data for instances
  for (odb::dbInst *db_inst : db_instances) {
    Instance instance(db_inst);
    instance.setDataStorage(&data_storage_);
    instance.setDataMapping(&data_mapping_);
    data_storage_.instances.push_back(instance);
  }
  // 2-3. make pointer set and map from db_instance to instance pointer.
  // Additionally: set the cell id
  for (int i = 0; i < data_storage_.instances.size(); ++i) {
    Instance *instance = &data_storage_.instances.at(i);
    instance_pointers_.push_back(instance);
    data_mapping_.inst_map[instance->getDbInst()] = instance;
    instance->setId(i);
  }


  /*!
   * @brief
   * Pin setting
   *
   * @details
   * Same with above way
   */
  // 1. make real data
  // 1-1. Instance terminals
  for (auto instance : instance_pointers_) {
    for (dbITerm *db_i_term : instance->getDbInst()->getITerms()) {
      Pin pin(db_i_term);
      pin.setDataStorage(&data_storage_);
      pin.setDataMapping(&data_mapping_);
      data_storage_.pins.push_back(pin);
    }
  }
  // 1-2. Block terminals
  for (dbBTerm *db_b_term : block->getBTerms()) {
    Pin pin(db_b_term);
    pin.setDataStorage(&data_storage_);
    pin.setDataMapping(&data_mapping_);
    data_storage_.pins.push_back(pin);
  }

  // 2. make pointer set and map from db_pin to pin pointer
  for (auto &pin : data_storage_.pins) {
    Pin *pin_pointer = &pin;
    pin_pointers_.push_back(pin_pointer);
    if (pin_pointer->isInstancePin()) {
      data_mapping_.pin_map_i[pin_pointer->getDbITerm()] = pin_pointer;
    } else if (pin_pointer->isBlockPin()) {
      data_mapping_.pin_map_b[pin_pointer->getDbBTerm()] = pin_pointer;
      pad_pointers_.push_back(pin_pointer);
    }
  }

  /*!
   * @brief
   * Net setting
   *
   * @details
   * Same with above way
   */
  data_storage_.nets.reserve(db_nets.size());
  net_pointers_.reserve(db_nets.size());
  // 1. make real data
  for (odb::dbNet *db_net : db_nets) {
    Net net(db_net);
    net.setDataStorage(&data_storage_);
    net.setDataMapping(&data_mapping_);
    data_storage_.nets.push_back(net);
  }
  // 2. make pointer set and map from db_net to net pointer
  for (auto &net : data_storage_.nets) {
    net_pointers_.push_back(&net);
    data_mapping_.net_map[net.getDbNet()] = &net;
  }


  /// Die setting
  // TODO: check whether this is valid
  int num_of_die = 2;
  data_storage_.dies.reserve(num_of_die);
  for (int i = 0; i < num_of_die; ++i) {
    Die die;
    die.setDbBlock(block);
    data_storage_.dies.push_back(die);
  }
  for (int i = 0; i < num_of_die; ++i) {
    die_pointers_.push_back(&data_storage_.dies.at(i));
  }

}
void Chip::write(const string &out_file_name) {
  parser_.writeDef(out_file_name);
}
ulong Chip::getHPWL() {
  ulong HPWL = 0;
  for (Net *net : net_pointers_) {
    HPWL += net->getHPWL();
  }
  return HPWL;
}
void Chip::parse_iccad(const string &lef_name, const string &def_name) {
/*
  // open input file
  ifstream input_file(input_file_name);
  if (input_file.fail()) {
    cerr << "Cannot open the input file: " << input_file_name << endl;
    exit(1);
  }

  // parsing start.
  string info, name1, name2;
  int n1, n2, n3, n4, n5;

  // check the number of Technologies
  // Syntax of input file: NumTechnologies <technologyCount>
  input_file >> info >> n1;
  assert(info == "NumTechnologies");
  this->num_technologies_ = n1;
  circuits_.resize(this->num_technologies_);

  dbTech *db_tech;
  dbTechLayer *db_tech_layer;
  dbLib *lib;
  dbChip *chip;
  dbBlock *block;

  // read Technologies
  vector<string> lib_cell_names;
  /////////////////////////////////// lef parsing for each tier ///////////////////////////////////////
  for (int i = 0; i < this->num_technologies_; ++i) {
    // Syntax: Tech <techName> <libCellCount>
    input_file >> info >> name1 >> n1;
    assert(info == "Tech");

    Chip *circuit = &circuits_.at(i);
    dbDatabase *db_database = circuit->getDbDatabase();

    // create tech
    // (refer to `dbDatabase* createSimpleDB()` submodule/OpenDB/tests/cpp/helper.cpp)
    db_tech = dbTech::create(circuit->getDbDatabase());
    db_tech_layer = dbTechLayer::create(db_tech, "metal1", dbTechLayerType::MASTERSLICE);
    lib = dbLib::create(db_database, info.c_str(), ',');
    chip = dbChip::create(db_database);
    block = dbBlock::create(chip, "block");

    // create libCells
    circuit->setLibCellNum(n1);
    for (int j = 0; j < circuit->getLibCellNum(); ++j) {
      // Syntax: LibCell <libCellName> <libCellSizeX> <libCellSizeY> <pinCount>
      input_file >> info >> name1 >> n2 >> n3 >> n1;
      assert(info == "LibCell");
      /// create one libCell
      // (refer to `dbDatabase* createSimpleDB()` submodule/OpenDB/tests/cpp/helper.cpp)
      dbMaster *master = dbMaster::create(lib, name1.c_str());
      master->setWidth(n2);
      master->setHeight(n3);
      master->setType(dbMasterType::CORE);
      if (i == 0)
        lib_cell_names.push_back(name1);
      /// create pins in the above libCell
      // (refer to addGeoms in `OpenDB/src/lefin/lefin.cpp`)
      for (int k = 0; k < n1; ++k) {
        // Syntax: Pin <pinName> <pinLocationX> <pinLocationY>
        // pin object generate
        input_file >> info >> name1 >> n2 >> n3;
        assert(info == "Pin");
        dbMTerm *db_m_term = dbMTerm::create(master, name1.c_str(), dbIoType::INOUT, dbSigType::SIGNAL);
        // pin geometry generate
        dbMPin *db_m_pin = dbMPin::create(db_m_term);
        int x1 = n2;
        int x2 = x1 + 1;
        int y1 = n3;
        int y2 = y1 + 1;
        dbBox::create(db_m_pin, db_tech_layer, x1, y2, x2, y2);
        // trim the order of pins in the master terminal
        dbSet<dbMPin> pins = db_m_term->getMPins();
        if (pins.reversible() && pins.orderReversed())
          pins.reverse();
      }
    }
  }

  // create pseudo lib cells for netlist
  db_tech = dbTech::create(db_database_netlist_);
  db_tech_layer = dbTechLayer::create(db_tech, "metal1", dbTechLayerType::MASTERSLICE);
  lib = dbLib::create(db_database_netlist_, "pseudo_layer", ',');
  chip = dbChip::create(db_database_netlist_);
  block = dbBlock::create(chip, "pseudo_block");
  for (int i = 0; i < lib_cell_names.size(); ++i) {
    dbMaster *master = dbMaster::create(lib, lib_cell_names.at(i).c_str());
    master->setType(dbMasterType::CORE);
  }

  // Syntax: DieSize <lowerLeftX> <lowerLeftY> <upperRightX> <upperRightY>
  input_file >> info >> n1 >> n2 >> n3 >> n4;
  assert(info == "DieSize");
  assert(n1 == 0 && n3 == 0);
  for (int i = 0; i < 2; ++i) {
    circuits_.at(i).setDieSize(static_cast<uint>(n2), static_cast<uint>(n4));
  }

  // Syntax: TopDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "TopDieMaxUtil");
  circuits_.at(0).setUtil(n1);

  // Syntax: BottomDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "BottomDieMaxUtil");
  circuits_.at(1).setUtil(n1);

  // Syntax: TopDieRows <startX> <startY> <rowLength> <rowHeight> <repeatCount>
  input_file >> info >> n2 >> n3 >> n4 >> n5 >> n1;
  assert(info == "TopDieRows");
  // TODO
  //  circuit.topDie.setRows(n2, n3, n4, n5, n1);

  // Syntax: BottomDieRows <startX> <startY> <rowLength> <rowHeight>
  // <repeatCount>
  input_file >> info >> n2 >> n3 >> n4 >> n5 >> n1;
  assert(info == "BottomDieRows");
  // TODO
  //  circuit.bottomDie.setRows(n2, n3, n4, n5, n1);

  // Syntax: TopDieTech <TechName>
  input_file >> info >> name1;
  assert(info == "TopDieTech");
  // TODO
  //  circuit.topDie.tech = techs[name1];
  //  for (auto &libCell : circuit.topDie.tech.libCells) {
  //    libCell.dieID = DieID::TOP;
  //  }

  // Syntax: BottomDieTech <TechName>
  input_file >> info >> name1;
  assert(info == "BottomDieTech");
  // TODO
  //  circuit.bottomDie.tech = techs[name1];
  //  for (auto &libCell : circuit.bottomDie.tech.libCells) {
  //    libCell.dieID = DieID::BOTTOM;
  //  }

  // Syntax: TerminalSize <sizeX> <sizeY>
  input_file >> info >> n1 >> n2;
  assert(info == "TerminalSize");
  this->hybrid_bond_size_.first = n1;
  this->hybrid_bond_size_.second = n2;

  // Syntax: TerminalSpacing <spacing>
  input_file >> info >> n1;
  assert(info == "TerminalSpacing");
  this->hybrid_bond_size_space_ = n1;

  // Syntax: NumInstances <instanceCount>
  input_file >> info >> n1;
  assert(info == "NumInstances");
  instance_num_ = n1;

  // read Instances in one circuit
  for (int i = 0; i < instance_num_; i++) {
    // Syntax: Inst <instName> <libCellName>
    input_file >> info >> name1 >> name2;
    assert(info == "Inst");
    auto master = lib->findMaster(name2.c_str());
    dbInst::create(block, master, name1.c_str());
  }

  // Syntax: NumNets <netCount>
  input_file >> info >> n1;
  assert(info == "NumNets");
  net_num_ = n1;

  // read Nets in one circuit
  // (refer to `create2LevetDbNoBTerms` function in submodule/OpenDB/tests/cpp/helper.cpp)
  for (int i = 0; i < net_num_; i++) {
    // Syntax: Net <netName> <numPins>
    input_file >> info >> name1 >> n1;
    assert(info == "Net");
    // TODO
    //    circuit.netlist.addNet(name1, n1);

    // read pins in one Net
    for (int j = 0; j < n1; j++) {
      // Syntax: Pin <instName>/<libPinName>
      input_file >> info >> name2;
      assert(info == "Pin");
      // TODO
      //      circuit.netlist.addPin(name1, name2);
    }
  }
*/
}
int Chip::getUnitOfMicro() const {
  return parser_.db_database_->getTech()->getDbUnitsPerMicron();
}

} // VLSI_backend