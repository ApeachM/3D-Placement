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
#include <fstream>
#include "Chip.h"

using namespace std;

namespace VLSI_backend {
Chip::Chip() {
}
void Chip::parse(const string &lef_name, const string &def_name) {
  parser_.readLef(lef_name);
  parser_.readDef(def_name);
  this->init();
}
void Chip::init() {
  /// this is for connection between the objects in initialization
  struct data_mapping {
    std::unordered_map<dbInst *, Instance *> inst_map;
    std::unordered_map<dbNet *, Net *> net_map;
    /// mapping for terminals on instance (pins on cell)
    std::unordered_map<dbITerm *, Pin *> pin_map_i;
    /// mapping for terminals on blocks (includes fixed pins on die)
    std::unordered_map<dbBTerm *, Pin *> pin_map_b;
  };

  dbBlock *block = db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();
  dbSet<dbNet> db_nets = block->getNets();
  data_mapping mapping;

  /**
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
  int id = 0;
  for (odb::dbInst *db_inst : db_instances) {
    Instance instance(db_inst, id);
    data_storage_.instances.push_back(instance);
    id++;
  }
  // 2-3. make pointer set into instance pointer.
  for (int i = 0; i < data_storage_.instances.size(); ++i) {
    Instance *instance = &data_storage_.instances.at(i);
    instance_pointers_.push_back(instance);
    mapping.inst_map[instance->getDbInst()] = instance;
  }

  /**
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
      data_storage_.pins.push_back(pin);
    }
  }
  // 1-2. Block terminals
  for (dbBTerm *db_b_term : block->getBTerms()) {
    Pin pin(db_b_term);
    data_storage_.pins.push_back(pin);
  }

  // 2. make pointer set into pin_pointers
  for (auto &pin : data_storage_.pins) {
    Pin *pin_pointer = &pin;
    pin_pointers_.push_back(pin_pointer);
    if (pin_pointer->isInstancePin()) {
      mapping.pin_map_i[pin_pointer->getDbITerm()] = pin_pointer;
    } else if (pin_pointer->isBlockPin()) {
      mapping.pin_map_b[pin_pointer->getDbBTerm()] = pin_pointer;
      pad_pointers_.push_back(pin_pointer);
    }
  }

  /**
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
    data_storage_.nets.push_back(net);
  }
  // 2. make pointer set and map from db_net to net pointer
  for (auto &net : data_storage_.nets) {
    Net *net_pointer = &net;
    net_pointers_.push_back(net_pointer);
    mapping.net_map[net_pointer->getDbNet()] = net_pointer;
  }

  /// Die setting
  // TODO: check whether this is valid
  int num_of_die = 2;
  data_storage_.dies.reserve(num_of_die);
  for (int i = 0; i < num_of_die + 1; ++i) {
    Die die;
    die.setDbBlock(block);
    die.setDieId(i);
    data_storage_.dies.push_back(die);
  }
  for (int i = 0; i < num_of_die + 1; ++i) {
    die_pointers_.push_back(&data_storage_.dies.at(i));
  }

  /// connect between [instance - pin - net]
  // instance -> pins  &&  instance -> nets
  for (Instance *instance : instance_pointers_) {
    vector<Net *> nets;
    vector<Pin *> pins;
    for (dbITerm *db_i_term : instance->getDbInst()->getITerms()) {
      Net *net = mapping.net_map[db_i_term->getNet()];
      Pin *pin = mapping.pin_map_i[db_i_term];
      nets.push_back(net);
      pins.push_back(pin);
    }
    instance->setConnectedNets(nets);
    instance->setConnectedPins(pins);
  }

  // pin -> instance  &&  pin -> net
  for (Pin *pin : pin_pointers_) {
    Instance *connected_instance;
    Net *connected_net;
    if (pin->isBlockPin()) {
      connected_instance = nullptr;
      connected_net = mapping.net_map[pin->getDbBTerm()->getNet()];
    } else if (pin->isInstancePin()) {
      connected_instance = mapping.inst_map[pin->getDbITerm()->getInst()];
      connected_net = mapping.net_map[pin->getDbITerm()->getNet()];
    }
    pin->setConnectedInstance(connected_instance);
    pin->setConnectedNet(connected_net);
  }

  // net -> instances  &&  net -> pins
  for (Net *net : net_pointers_) {
    vector<Instance *> instances;
    vector<Pin *> pins;
    for (dbITerm *i_term : net->getDbNet()->getITerms()) {
      Instance *instance = mapping.inst_map[i_term->getInst()];
      Pin *pin = mapping.pin_map_i[i_term];
      instances.push_back(instance);
      pins.push_back(pin);
    }
    for (dbBTerm *b_term : net->getDbNet()->getBTerms()) {
      Pin *pin = mapping.pin_map_b[b_term];
      pins.push_back(pin);
    }
    net->setConnectedInstances(instances);
    net->setConnectedPins(pins);
  }

  /* Warning! Never use `push_back` method for `data_storage_.[something]` after this line.
  *
  * Recommend: If you want to call `push_back` for `data_storage_.[something]` inevitably,
  * then just make another data_storage for objects and
  * push back their pointers into `~_pointers` (i.e. instance_pointers).
  * Then there will be no problem.
  * */

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
void Chip::parseICCAD(const string &input_file_name) {
  // In this function, we only construct odb database.

  // open input file
  ifstream input_file(input_file_name);
  if (input_file.fail()) {
    cerr << "Cannot open the input file: " << input_file_name << endl;
    exit(1);
  }

  // These are for collecting top and bottom ones and making pseudo die library with these infos
  struct PinInfo {
    string pin_name;
    int pin_location_x;
    int pin_location_y;
  };
  struct LibCellInfo {
    string name;
    int width;
    int height;
    int pin_number;
    vector<PinInfo> pin_infos;
  };
  vector<LibCellInfo> lib_cell_infos_top;
  vector<LibCellInfo> lib_cell_infos_bottom;


  // parsing start //

  // temporal variables
  string info, name1, name2;
  int n1, n2, n3, n4, n5;

  // check the number of Technologies
  // Syntax of input file: NumTechnologies <technologyCount>
  input_file >> info >> n1;
  assert(info == "NumTechnologies");
  num_technologies_ = n1;

  // Library parsing(lef parsing) for each tier
  for (int die_id = 0; die_id < num_technologies_; ++die_id) {
    // Syntax: Tech <techName> <libCellCount>
    input_file >> info >> name1 >> n1;
    assert(info == "Tech");

    assert(db_databases_.size() == die_id);
    dbDatabase *db_database = dbDatabase::create();
    dbTech *db_tech = dbTech::create(db_database);
    dbTechLayer *db_tech_layer = dbTechLayer::create(db_tech, "Layer", dbTechLayerType::MASTERSLICE);
    dbLib *db_lib = dbLib::create(db_database, name1.c_str(), ',');
    dbChip *db_chip = dbChip::create(db_database);
    dbBlock *db_block = dbBlock::create(db_chip, (std::to_string(die_id) + "th Die Block").c_str());
    db_databases_.push_back(db_database);

    // read LibCells in one tech
    for (int i = 0; i < n1; ++i) {
      // Syntax: LibCell <libCellName> <libCellSizeX> <libCellSizeY> <pinCount>
      input_file >> info >> name1 >> n1 >> n2 >> n3;
      assert(info == "LibCell");

      LibCellInfo lib_cell_info;
      lib_cell_info.name = name1;
      lib_cell_info.width = n1;
      lib_cell_info.height = n2;
      lib_cell_info.pin_number = n3;

      // (refer to `dbDatabase* createMaster2X1()` submodule/OpenDB/tests/cpp/helper.cpp)
      dbMaster *master = dbMaster::create(db_lib, name1.c_str());
      master->setWidth(lib_cell_info.width);
      master->setHeight(lib_cell_info.width);
      master->setType(dbMasterType::CORE);
      // read pins in one LibCell
      for (int j = 0; j < lib_cell_info.pin_number; ++j) {
        // Syntax: Pin <pinName> <pinLocationX> <pinLocationY>
        input_file >> info >> name1 >> n4 >> n5;
        assert(info == "Pin");
        PinInfo pin_info;
        pin_info.pin_name = name1;
        pin_info.pin_location_x = n4;
        pin_info.pin_location_y = n5;
        lib_cell_info.pin_infos.push_back(pin_info);

        // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
        dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
        dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
        dbMTerm *master_terminal = dbMTerm::create(master, name1.c_str(), io_type, sig_type);
        dbMPin *db_m_pin = dbMPin::create(master_terminal);
        // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
        dbBox::create(db_m_pin, db_tech_layer,
                      pin_info.pin_location_x,
                      pin_info.pin_location_y,
                      pin_info.pin_location_x + 1,
                      pin_info.pin_location_y + 1);
      }
      if (die_id == 0)
        lib_cell_infos_top.push_back(lib_cell_info);
      else if (die_id == 1)
        lib_cell_infos_bottom.push_back(lib_cell_info);
    }
  }


  // Library parsing(lef parsing) for pseudo tier
  assert(db_database_ == nullptr);
  db_database_ = odb::dbDatabase::create();
  dbTech *db_tech = dbTech::create(db_database_);
  dbTechLayer *db_tech_layer = dbTechLayer::create(db_tech, "pseudoLayer", dbTechLayerType::MASTERSLICE);
  dbLib *db_lib = dbLib::create(db_database_, "pseudoDieLib", ',');
  dbChip *db_chip = dbChip::create(db_database_);
  dbBlock *db_block = dbBlock::create(db_chip, "Pseudo Die Block");

  assert(lib_cell_infos_top.size() == lib_cell_infos_bottom.size());
  for (int i = 0; i < lib_cell_infos_top.size(); ++i) {
    string lib_cell_name;
    int width, height;
    LibCellInfo lib_cell_info_top = lib_cell_infos_top.at(i);
    LibCellInfo lib_cell_info_bottom = lib_cell_infos_bottom.at(i);
    assert(lib_cell_info_top.name == lib_cell_info_bottom.name);
    lib_cell_name = lib_cell_info_top.name;
    width = floor((lib_cell_info_top.width + lib_cell_info_bottom.width) / 2);
    height = floor((lib_cell_info_top.height + lib_cell_info_bottom.height) / 2);

    dbMaster *master = dbMaster::create(db_lib, lib_cell_name.c_str());
    master->setWidth(width);
    master->setHeight(height);
    master->setType(dbMasterType::CORE);

    assert(lib_cell_info_top.pin_number == lib_cell_info_bottom.pin_number);
    for (int j = 0; j < lib_cell_info_top.pin_number; ++j) {
      PinInfo pin_info_top = lib_cell_info_top.pin_infos.at(j);
      PinInfo pin_info_bottom = lib_cell_info_bottom.pin_infos.at(j);
      assert(pin_info_top.pin_name == pin_info_bottom.pin_name);
      string pin_name = pin_info_top.pin_name;
      int pin_location_x = floor((pin_info_top.pin_location_x + pin_info_bottom.pin_location_x) / 2);
      int pin_location_y = floor((pin_info_top.pin_location_y + pin_info_bottom.pin_location_y) / 2);
      assert(width > pin_location_x);
      assert(height > pin_location_y);

      // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
      dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
      dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
      dbMTerm *master_terminal = dbMTerm::create(master, pin_name.c_str(), io_type, sig_type);
      dbMPin *db_m_pin = dbMPin::create(master_terminal);
      // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
      dbBox::create(db_m_pin, db_tech_layer, pin_location_x, pin_location_y, pin_location_x + 1, pin_location_y + 1);
    }
  }

  // Syntax: DieSize <lowerLeftX> <lowerLeftY> <upperRightX> <upperRightY>
  input_file >> info >> n1 >> n2 >> n3 >> n4;
  assert(info == "DieSize");
  odb::Point lower_left = odb::Point(n1, n2);
  odb::Point upper_right = odb::Point(n3, n4);
  odb::Rect rect(lower_left.getX(), lower_left.getY(), upper_right.getX(), upper_right.getY());

  // refer to submodule/OpenDB/src//defin/definReader.cpp
  db_database_->getChip()->getBlock()->setDieArea(rect);
  db_databases_.at(0)->getChip()->getBlock()->setDieArea(rect);
  db_databases_.at(0)->getChip()->getBlock()->setDieArea(rect);

  // Syntax: TopDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "TopDieMaxUtil");
  max_util_.first = n1;

  // Syntax: BottomDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "BottomDieMaxUtil");
  max_util_.second = n1;

  // Syntax: TerminalSize <sizeX> <sizeY>
  input_file >> info >> n1 >> n2;
  assert(info == "TerminalSize");
  hybrid_size_x_ = n1;
  hybrid_size_y_ = n2;

  // Syntax: TerminalSpacing <spacing>
  input_file >> info >> n1;
  assert(info == "TerminalSpacing");
  hybrid_spacing_ = n1;

  // Syntax: NumInstances <instanceCount>
  input_file >> info >> n1;
  assert(info == "NumInstances");
  instance_number_ = n1;

  // read Instances in one circuit
  for (int i = 0; i < instance_number_; ++i) {
    // Syntax: Inst <instName> <libCellName>
    input_file >> info >> name1 >> name2;
    assert(info == "Inst");
    dbMaster *master = db_database_->findMaster(name2.c_str());
    dbInst::create(db_block, master, name1.c_str());
  }

  // Syntax: NumNets <netCount>
  input_file >> info >> n1;
  assert(info == "NumNets");
  net_number_ = n1;

  // read Nets in one circuit
  for (int i = 0; i < net_number_; ++i) {
    // (refer to `dbDatabase* create2LevetDbNoBTerms()` function in submodule/OpenDB/test/cpp/helper.cpp)
    // Syntax: Net <netName> <numPins>
    input_file >> info >> name1 >> n1;
    assert(info == "Net");
    dbNet* net = dbNet::create(db_block, name1.c_str());

    // read pins in one Net
    for (int j = 0; j < n1; ++j) {
      // Syntax: Pin <instName>/<libPinName>
      input_file >> info >> name2;
      assert(info == "Pin");

      int idx = name2.find('/');
      string inst_name =  name2.substr(0, idx);
      string lib_pin_name = name2.substr(idx+1);

      dbInst* inst = db_block->findInst(inst_name.c_str());
      dbITerm::connect(inst->findITerm(lib_pin_name.c_str()), net);
    }
  }


  this->init();
}
void Chip::parseICCAD_deprecated(const string &input_file_name) {
  struct PinInfo {
    string pin_name;
    int pin_location_x;
    int pin_location_y;
  };
  struct LibCellInfo {
    string lib_cell_name;
    int lib_cell_size_x;
    int lib_cell_size_y;
    int pin_number;
    vector<PinInfo> pin_infos;
  };

  vector<LibCellInfo> lib_cell_infos1;
  vector<LibCellInfo> lib_cell_infos2;

  // open input file
  ifstream input_file(input_file_name);
  if (input_file.fail()) {
    cerr << "Cannot open the input file: " << input_file_name << endl;
    exit(1);
  }
  // parsing start //

  // temporal variables
  string info, name1, name2;
  int n1, n2, n3, n4, n5;

  // check the number of Technologies
  // Syntax of input file: NumTechnologies <technologyCount>
  input_file >> info >> n1;
  assert(info == "NumTechnologies");
  num_technologies_ = n1;

  // Make Dies as much as technologies
  // one for pseudo die
  for (int i = 0; i < num_technologies_ + 1; ++i) {
    Die die;
    die.setDieId(i);
    data_storage_.dies.push_back(die);
  }
  for (int i = 0; i < n1 + 1; ++i)
    die_pointers_.push_back(&data_storage_.dies.at(i));

  /////////////////////////////////// lef parsing for each tier ///////////////////////////////////////
  for (int i = 1; i < num_technologies_ + 1; ++i) {
    // Syntax: Tech <techName> <libCellCount>
    input_file >> info >> name1 >> n1;
    assert(info == "Tech");

    // Tech tech(name1, n1);
    Die *target_die = die_pointers_.at(i);
    target_die->setDBBasic(name1);
    dbDatabase *db_database = target_die->getDbDatabase();
    dbTech *db_tech = target_die->getDbTech();
    dbTechLayer *db_tech_layer = target_die->getDbTechLayer();
    dbLib *db_lib = target_die->getDbLib();
    dbChip *db_chip = target_die->getDbChip();
    dbBlock *db_block = target_die->getDbBlock();
    target_die->setLibNum(n1);

    // read LibCells in one tech
    for (int j = 0; j < die_pointers_.at(i)->getLibNum(); j++) {
      // Syntax: LibCell <libCellName> <libCellSizeX> <libCellSizeY> <pinCount>
      input_file >> info >> name1 >> n1 >> n2 >> n3;
      assert(info == "LibCell");

      LibCellInfo lib_cell_info;
      lib_cell_info.lib_cell_name = name1;
      lib_cell_info.lib_cell_size_x = n1;
      lib_cell_info.lib_cell_size_y = n2;
      lib_cell_info.pin_number = n3;

      // (refer to `dbDatabase* createMaster2X1()` submodule/OpenDB/tests/cpp/helper.cpp)
      dbMaster *master = dbMaster::create(db_lib, name1.c_str());
      master->setWidth(n1);
      master->setHeight(n2);
      master->setType(dbMasterType::CORE);
      // read pins in one LibCell
      for (int k = 0; k < n3; ++k) {
        // Syntax: Pin <pinName> <pinLocationX> <pinLocationY>
        input_file >> info >> name1 >> n4 >> n5;
        assert(info == "Pin");

        PinInfo pin_info;
        pin_info.pin_name = name1;
        pin_info.pin_location_x = n4;
        pin_info.pin_location_y = n5;
        lib_cell_info.pin_infos.push_back(pin_info);

        // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
        dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
        dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
        dbMTerm *master_terminal = dbMTerm::create(master, name1.c_str(), io_type, sig_type);
        dbMPin *db_m_pin = dbMPin::create(master_terminal);

        // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
        dbBox::create(db_m_pin, db_tech_layer, n4, n5, n4 + 1, n5 + 1);
      }
      if (i == 1)
        lib_cell_infos1.push_back(lib_cell_info);
      else if (i == 2)
        lib_cell_infos2.push_back(lib_cell_info);
    }
  }

  /////////////////////////////////// lef parsing for pseudo tier ///////////////////////////////////////
  die_pointers_.at(0)->setDbDatabase(odb::dbDatabase::create());
  dbDatabase *db_database = die_pointers_.at(0)->getDbDatabase();
  dbTech *db_tech = dbTech::create(db_database);
  dbTechLayer *db_tech_layer = dbTechLayer::create(db_tech, "pseudoLayer", dbTechLayerType::MASTERSLICE);
  dbLib *db_lib = dbLib::create(db_database, "pseudoDieLib", ',');
  dbChip *db_chip = dbChip::create(db_database);
  dbBlock *db_block = dbBlock::create(db_chip, (to_string(0) + "th Die Block").c_str());

  die_pointers_.at(0)->setDbTech(db_tech);
  die_pointers_.at(0)->setDbTechLayer(db_tech_layer);
  die_pointers_.at(0)->setDbLib(db_lib);
  die_pointers_.at(0)->setDbChip(db_chip);
  die_pointers_.at(0)->setDbBlock(db_block);

  // set LibCells for pseudo die
  for (int i = 0; i < lib_cell_infos1.size(); ++i) {
    string lib_cell_name;
    int width, height;
    LibCellInfo lib_cell_info1 = lib_cell_infos1.at(i);
    LibCellInfo lib_cell_info2 = lib_cell_infos2.at(i);
    assert(lib_cell_info1.lib_cell_name == lib_cell_info2.lib_cell_name);
    lib_cell_name = lib_cell_info1.lib_cell_name;
    width = floor((lib_cell_info1.lib_cell_size_x + lib_cell_info2.lib_cell_size_x) / 2);
    height = floor((lib_cell_info1.lib_cell_size_y + lib_cell_info2.lib_cell_size_y) / 2);

    dbMaster *master = dbMaster::create(db_lib, lib_cell_name.c_str());
    master->setWidth(width);
    master->setHeight(height);
    master->setType(dbMasterType::CORE);

    assert(lib_cell_info1.pin_number == lib_cell_info2.pin_number);
    for (int j = 0; j < lib_cell_info1.pin_number; ++j) {
      PinInfo pin_info1 = lib_cell_info1.pin_infos.at(j);
      PinInfo pin_info2 = lib_cell_info2.pin_infos.at(j);
      assert(pin_info1.pin_name == pin_info2.pin_name);
      string pin_name = pin_info1.pin_name;
      int pin_location_x = floor((pin_info1.pin_location_x + pin_info2.pin_location_x) / 2);
      int pin_location_y = floor((pin_info1.pin_location_y + pin_info2.pin_location_y) / 2);
      assert(width > pin_location_x);
      assert(height > pin_location_y);

      // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
      dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
      dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
      dbMTerm *master_terminal = dbMTerm::create(master, pin_name.c_str(), io_type, sig_type);
      dbMPin *db_m_pin = dbMPin::create(master_terminal);

      // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
      dbBox::create(db_m_pin, db_tech_layer, pin_location_x, pin_location_y, pin_location_x + 1, pin_location_y + 1);
    }
  }

  // Syntax: DieSize <lowerLeftX> <lowerLeftY> <upperRightX> <upperRightY>
  input_file >> info >> n1 >> n2 >> n3 >> n4;
  assert(info == "DieSize");
  odb::Point lower_left = odb::Point(n1, n2);
  odb::Point upper_right = odb::Point(n3, n4);
  odb::Rect rect(lower_left.getX(), lower_left.getY(), upper_right.getX(), upper_right.getY());
  for (int i = 0; i < num_technologies_ + 1; ++i) {
    // refer to submodule/OpenDB/src//defin/definReader.cpp
    die_pointers_.at(i)->getDbBlock()->setDieArea(rect);
  }

  // Syntax: TopDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "TopDieMaxUtil");
  die_pointers_.at(1)->setMaxUtil(n1);

  // Syntax: BottomDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "BottomDieMaxUtil");
  die_pointers_.at(2)->setMaxUtil(n1);

  // for pseudo die
  die_pointers_.at(0)->setMaxUtil(100);

  // Syntax: TerminalSize <sizeX> <sizeY>
  input_file >> info >> n1 >> n2;
  assert(info == "TerminalSize");
  hybrid_size_x_ = n1;
  hybrid_size_y_ = n2;

  // Syntax: TerminalSpacing <spacing>
  input_file >> info >> n1;
  assert(info == "TerminalSpacing");
  hybrid_spacing_ = n1;

  // Syntax: NumInstances <instanceCount>
  input_file >> info >> n1;
  assert(info == "NumInstances");
  instance_number_ = n1;

  // read Instances in one circuit
  for (int i = 0; i < instance_number_; ++i) {
    // Syntax: Inst <instName> <libCellName>
    input_file >> info >> name1 >> name2;
    assert(info == "Inst");
    Instance instance;
    instance.setInstName(name1);
    instance.setLibName(name2);
    dbMaster *master = db_database->findMaster(instance.getLibName().c_str());
    instance.setLibrary(master);
  }

  // read Instances in pseudo circuit
  // The other 2 die will be used after finishing all process and before writing
  for (int i = 0; i < instance_number_; ++i) {
    // Syntax: Inst <instName> <libCellName>
    input_file >> info >> name1 >> name2;
    assert(info == "Inst");
    Instance instance;
    instance.setInstName(name1);
    instance.setLibName(name2);
    dbMaster *master = db_lib->findMaster(instance.getLibName().c_str());
    instance.setLibrary(master);
    data_storage_.instances.push_back(instance);
  }
  // instance pointer setting
  assert(instance_pointers_.empty());
  for (int i = 0; i < instance_number_; ++i) {
    Instance *instance = &data_storage_.instances.at(i);
    instance_pointers_.push_back(instance);
  }

  // Syntax: NumNets <netCount>
  input_file >> info >> n1;
  assert(info == "NumNets");
  net_number_ = n1;

  // read Nets in one circuit
  for (int i = 0; i < net_number_; ++i) {
    // Syntax: Net <netName> <numPins>
    input_file >> info >> name1 >> n1;
    assert(info == "Net");

    // read pins in one Net
    for (int j = 0; j < n1; ++j) {
      // Syntax: Pin <instName>/<libPinName>
      input_file >> info >> name2;
      assert(info == "Pin");

    }
  }

}
void Chip::test() {
  dbDatabase *die_db = odb::dbDatabase::create();
  setDbDatabase(die_db);

  // create simple db
  cout << die_db->getTech() << endl;

  dbTech *db_tech = dbTech::create(die_db);
  dbTechLayer *layer = dbTechLayer::create(db_tech, "pseudoLayer1", dbTechLayerType::MASTERSLICE);
  dbLib *lib = dbLib::create(die_db, "lib", ',');
  dbChip *chip = dbChip::create(die_db);
  dbBlock *block = dbBlock::create(chip, "simple_block");

  dbMaster *master_and2 = dbMaster::create(lib, "and2");
  master_and2->setWidth(30);
  master_and2->setHeight(40);
  master_and2->setType(dbMasterType::CORE);
  dbMTerm *master_terminal_a = dbMTerm::create(master_and2, "a", dbIoType::INOUT, dbSigType::SIGNAL);
  dbMTerm *master_terminal_b = dbMTerm::create(master_and2, "b", dbIoType::INPUT, dbSigType::SIGNAL);
  dbMTerm *master_terminal_o = dbMTerm::create(master_and2, "o", dbIoType::OUTPUT, dbSigType::SIGNAL);
  dbMPin *m_pin_a = dbMPin::create(master_terminal_a);
  dbMPin *m_pin_b = dbMPin::create(master_terminal_b);
  dbMPin *m_pin_o = dbMPin::create(master_terminal_o);
  dbBox::create(m_pin_a, layer, 1, 2, 3, 4);
  dbBox::create(m_pin_b, layer, 10, 20, 11, 21);
  dbBox::create(m_pin_o, layer, 20, 30, 21, 31);
  master_and2->setFrozen();

  dbMaster *master_or2 = dbMaster::create(lib, "or");
  master_or2->setWidth(40);
  master_or2->setHeight(50);
  master_or2->setType(dbMasterType::CORE);
  master_terminal_a = dbMTerm::create(master_or2, "a", dbIoType::INOUT, dbSigType::SIGNAL);
  master_terminal_b = dbMTerm::create(master_or2, "b", dbIoType::INPUT, dbSigType::SIGNAL);
  master_terminal_o = dbMTerm::create(master_or2, "o", dbIoType::OUTPUT, dbSigType::SIGNAL);
  m_pin_a = dbMPin::create(master_terminal_a);
  m_pin_b = dbMPin::create(master_terminal_b);
  m_pin_o = dbMPin::create(master_terminal_o);
  dbBox::create(m_pin_a, layer, 2, 3, 3, 4);
  dbBox::create(m_pin_b, layer, 10, 20, 11, 21);
  dbBox::create(m_pin_o, layer, 30, 40, 31, 41);
  master_or2->setFrozen();


  // #     (n1)   +-----
  // #    --------|a    \    (n5)
  // #     (n2)   | (i1)o|-----------+
  // #    --------|b    /            |       +-------
  // #            +-----             +--------\a     \    (n7)
  // #                                         ) (i3)o|---------------
  // #     (n3)   +-----             +--------/b     /
  // #    --------|a    \    (n6)    |       +-------
  // #     (n4)   | (i2)o|-----------+
  // #    --------|b    /
  // #            +-----

  // create 2 level db with no b terms
  dbInst *i1 = dbInst::create(block, master_and2, "i1");
  dbInst *i2 = dbInst::create(block, master_and2, "i2");
  dbInst *i3 = dbInst::create(block, master_or2, "i3");

  dbNet *n1 = dbNet::create(block, "n1");
  dbNet *n2 = dbNet::create(block, "n2");
  dbNet *n3 = dbNet::create(block, "n3");
  dbNet *n4 = dbNet::create(block, "n4");
  dbNet *n5 = dbNet::create(block, "n5");
  dbNet *n6 = dbNet::create(block, "n6");
  dbNet *n7 = dbNet::create(block, "n7");
  dbITerm::connect(i1->findITerm("a"), n1);
  dbITerm::connect(i1->findITerm("b"), n2);
  dbITerm::connect(i2->findITerm("a"), n3);
  dbITerm::connect(i2->findITerm("b"), n4);
  dbITerm::connect(i3->findITerm("a"), n5);
  dbITerm::connect(i3->findITerm("b"), n6);
  dbITerm::connect(i1->findITerm("o"), n5);
  dbITerm::connect(i2->findITerm("o"), n6);
  dbITerm::connect(i3->findITerm("o"), n7);

  i1->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i1->setLocation(30, 20);
  i2->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i2->setLocation(200, 300);
  i3->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i3->setLocation(0, 0);

  cout << "chip: " << chip << "\t" << db_database_->getChip() << endl;
  cout << "block: " << block << "\t" << db_database_->getChip()->getBlock() << endl;

}
int Chip::getUnitOfMicro() const {
  return db_database_->getTech()->getDbUnitsPerMicron();
}

} // VLSI_backend