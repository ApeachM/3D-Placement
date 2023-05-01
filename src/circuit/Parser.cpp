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

#include <fstream>
#include "Chip.h"

namespace VLSI_backend {
void Chip::parse(const string &lef_name, const string &def_name) {
  db_database_ = dbDatabase::create();
  parser_.db_database_ = db_database_;
  parser_.readLef(lef_name);
  parser_.readDef(def_name);
  this->init();
}
void Chip::write(const string &out_file_name) {
  parser_.writeDef(out_file_name);
}
void Chip::parseICCAD(const string &input_file_name) {
  struct LibPinInfo {
    string pin_name;
    int pin_location_x;
    int pin_location_y;
  };
  struct LibCellInfo {
    string name;
    int width;
    int height;
    int pin_number;
    vector<LibPinInfo> lib_pin_infos;
  };
  struct TechInfo {
    string name;
    int lib_cell_num;
    vector<LibCellInfo> lib_cell_infos;
  };
  struct DieInfo {
    string tech_name;
    int lower_left_x;
    int lower_left_y;
    int upper_right_x;
    int upper_right_y;
    int max_util;
    RowInfo row_info;
    TechInfo *tech_info;
  };
  struct TerminalInfo {
    int size_x;
    int size_y;
    int spacing_size;
  };
  struct InstanceInfo {
    string inst_name;
    string lib_cell_name;
  };
  struct ConnectedPinInfo {
    string instance_name;
    string lib_pin_name;
  };
  struct NetInfo {
    string net_name;
    int pin_num;
    vector<ConnectedPinInfo> connected_pins_infos;
  };

  vector<TechInfo> tech_infos;
  vector<InstanceInfo> instance_infos;
  vector<NetInfo> net_infos;
  vector<DieInfo> die_infos;
  TerminalInfo terminal_info{};


  // In this function, we only construct odb database.
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

  // Syntax: NumTechnologies <technologyCount>
  input_file >> info >> n1;
  assert(info == "NumTechnologies");
  num_technologies_ = n1;

  for (int i = 0; i < num_technologies_; ++i) {
    // Syntax: Tech <techName> <libCellCount>
    input_file >> info >> name1 >> n1;
    assert(info == "Tech");

    TechInfo tech_info;
    tech_info.name = name1;
    tech_info.lib_cell_num = n1;

    // read LibCells in one tech
    for (int j = 0; j < tech_info.lib_cell_num; ++j) {
      // Syntax: LibCell <libCellName> <libCellSizeX> <libCellSizeY> <pinCount>
      input_file >> info >> name1 >> n1 >> n2 >> n3;
      assert(info == "LibCell");

      LibCellInfo lib_cell_info;
      lib_cell_info.name = name1;
      lib_cell_info.width = n1;
      lib_cell_info.height = n2;
      lib_cell_info.pin_number = n3;

      for (int k = 0; k < lib_cell_info.pin_number; ++k) {
        // Syntax: Pin <pinName> <pinLocationX> <pinLocationY>
        input_file >> info >> name1 >> n1 >> n2;
        LibPinInfo lib_pin_info;
        lib_pin_info.pin_name = name1;
        lib_pin_info.pin_location_x = n1;
        lib_pin_info.pin_location_y = n2;
        lib_cell_info.lib_pin_infos.push_back(lib_pin_info);
      }
      tech_info.lib_cell_infos.push_back(lib_cell_info);
    }
    tech_infos.push_back(tech_info);
  }

  // Syntax: DieSize <lowerLeftX> <lowerLeftY> <upperRightX> <upperRightY>
  input_file >> info >> n1 >> n2 >> n3 >> n4;
  assert(info == "DieSize");

  for (int i = 0; i < 2; ++i) {
    DieInfo die_info;
    die_info.lower_left_x = n1;
    die_info.lower_left_y = n2;
    die_info.upper_right_x = n3;
    die_info.upper_right_y = n4;
    die_infos.push_back(die_info);
  }

  // Syntax: TopDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "TopDieMaxUtil");
  die_infos.at(0).max_util = n1;

  // Syntax: BottomDieMaxUtil <util>
  input_file >> info >> n1;
  assert(info == "BottomDieMaxUtil");
  die_infos.at(1).max_util = n1;

  // Syntax: TopDieRows <startX> <startY> <rowLength> <rowHeight> <repeatCount>
  input_file >> info >> n1 >> n2 >> n3 >> n4 >> n5;
  assert(info == "TopDieRows");
  die_infos.at(0).row_info.start_x = n1;
  die_infos.at(0).row_info.start_y = n2;
  die_infos.at(0).row_info.row_width = n3;
  die_infos.at(0).row_info.row_height = n4;
  die_infos.at(0).row_info.repeat_count = n5;

  // Syntax: BottomDieRows <startX> <startY> <rowLength> <rowHeight> <repeatCount>
  input_file >> info >> n1 >> n2 >> n3 >> n4 >> n5;
  assert(info == "BottomDieRows");
  die_infos.at(1).row_info.start_x = n1;
  die_infos.at(1).row_info.start_y = n2;
  die_infos.at(1).row_info.row_width = n3;
  die_infos.at(1).row_info.row_height = n4;
  die_infos.at(1).row_info.repeat_count = n5;

  // Syntax: TopDieTech <TechName>
  input_file >> info >> name1;
  assert(info == "TopDieTech");
  die_infos.at(0).tech_name = name1;

  // Syntax: BottomDieTech <TechName>
  input_file >> info >> name1;
  assert(info == "BottomDieTech");
  die_infos.at(1).tech_name = name1;

  // Syntax: TerminalSize <sizeX> <sizeY>
  input_file >> info >> n1 >> n2;
  assert(info == "TerminalSize");
  terminal_info.size_x = n1;
  terminal_info.size_y = n2;

  // Syntax: TerminalSpacing <spacing>
  input_file >> info >> n1;
  assert(info == "TerminalSpacing");
  terminal_info.spacing_size = n1;

  // Syntax: NumInstances <instanceCount>
  input_file >> info >> n1;
  assert(info == "NumInstances");
  instance_number_ = n1;

  // read instances
  for (int i = 0; i < instance_number_; ++i) {
    // Syntax: Inst <instName> <libCellName>
    input_file >> info >> name1 >> name2;
    assert(info == "Inst");
    InstanceInfo instance_info;
    instance_info.inst_name = name1;
    instance_info.lib_cell_name = name2;
    instance_infos.push_back(instance_info);
  }

  // Syntax: NumNets <netCount>
  input_file >> info >> n1;
  assert(info == "NumNets");
  net_number_ = n1;

  for (int i = 0; i < net_number_; ++i) {
    // Syntax: Net <netName> <numPins>
    input_file >> info >> name1 >> n1;
    assert(info == "Net");
    NetInfo net_info;
    net_info.net_name = name1;
    net_info.pin_num = n1;

    // read pins in one Net
    for (int j = 0; j < net_info.pin_num; ++j) {
      // Syntax: Pin <instName>/<libPinName>
      input_file >> info >> name1;
      assert(info == "Pin");

      int idx = name1.find('/');
      string inst_name = name1.substr(0, idx);
      string lib_pin_name = name1.substr(idx + 1);
      ConnectedPinInfo pin_info;
      pin_info.instance_name = inst_name;
      pin_info.lib_pin_name = lib_pin_name;
      net_info.connected_pins_infos.push_back(pin_info);
    }
    net_infos.push_back(net_info);
  }

  for (DieInfo &die_info : die_infos)
    for (TechInfo &tech_info : tech_infos)
      if (die_info.tech_name == tech_info.name)
        die_info.tech_info = &tech_info;

  row_infos_.first = die_infos.at(0).row_info;
  max_utils_.first = die_infos.at(0).max_util;
  row_infos_.second = die_infos.at(1).row_info;
  max_utils_.second = die_infos.at(1).max_util;

  // check
  for (DieInfo die_info : die_infos) {
    assert(die_info.tech_name == die_info.tech_info->name);
  }


  // Data parsing is completed.
  // Now, we construct odb database.

  // DB Base Construction //
  // for top and bottom die
  assert(db_databases_.empty());
  for (int i = 0; i < 2; ++i) {
    string die_name;
    if (i == 0)
      die_name = "Top Die";
    else if (i == 1)
      die_name = "Bottom Die";
    dbDatabase *db_database = odb::dbDatabase::create();
    dbTech *db_tech = dbTech::create(db_database);
    dbTechLayer::create(db_tech, die_name.c_str(), dbTechLayerType::MASTERSLICE);
    dbLib::create(db_database, "lib", ',');
    dbChip *db_chip = dbChip::create(db_database);
    dbBlock::create(db_chip, (die_name + " Block").c_str());
    db_databases_.push_back(db_database);
  }

  // for pseudo die
  {
    assert(db_database_ == nullptr);
    db_database_ = odb::dbDatabase::create();
    db_database_->setLogger(&logger_);
    dbTech *db_tech = dbTech::create(db_database_);
    dbTechLayer::create(db_tech, "pseudoLayer", dbTechLayerType::MASTERSLICE);
    dbLib::create(db_database_, "pseudoDieLib", ',');
    dbChip *db_chip = dbChip::create(db_database_);
    dbBlock::create(db_chip, "Pseudo Die Block");
  }


  // Library Construction //
  // for top and bottom die
  for (int die_id = 0; die_id < 2; ++die_id) {
    DieInfo *die_info = &die_infos.at(die_id);
    TechInfo *tech_info = die_info->tech_info;
    dbDatabase *db_database = db_databases_.at(die_id);
    dbTech *db_tech = db_database->getTech();
    dbTechLayer *db_tech_layer = db_tech->findLayer(0);
    dbLib *db_lib = db_database->findLib("lib");
    dbChip *db_chip = db_database->getChip();
    dbBlock *db_block = db_chip->getBlock();
    assert(db_tech_layer != nullptr);

    for (int i = 0; i < tech_info->lib_cell_num; ++i) {
      LibCellInfo lib_cell_info = tech_info->lib_cell_infos.at(i);
      // (refer to `dbDatabase* createMaster2X1()` submodule/OpenDB/tests/cpp/helper.cpp)
      dbMaster *master = dbMaster::create(db_lib, lib_cell_info.name.c_str());
      master->setWidth(lib_cell_info.width);
      master->setHeight(lib_cell_info.height);
      master->setType(dbMasterType::CORE);
      // read pins in one LIbCell
      for (int j = 0; j < lib_cell_info.pin_number; ++j) {
        LibPinInfo pin_info = lib_cell_info.lib_pin_infos.at(j);
        // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
        dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
        dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
        dbMTerm *master_terminal = dbMTerm::create(master, pin_info.pin_name.c_str());
        dbMPin *db_m_pin = dbMPin::create(master_terminal);
        // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
        dbBox::create(db_m_pin, db_tech_layer,
                      pin_info.pin_location_x,
                      pin_info.pin_location_y,
                      pin_info.pin_location_x + 1,
                      pin_info.pin_location_y + 1);
      }
      master->setFrozen();
    }
  }
  // for pseudo die
  dbTech *db_tech = db_database_->getTech();
  dbTechLayer *db_tech_layer = db_tech->findLayer(0);
  dbLib *db_lib = db_database_->findLib("pseudoDieLib");
  dbChip *db_chip = db_database_->getChip();
  dbBlock *db_block = db_chip->getBlock();
  assert(db_tech_layer != nullptr);

  DieInfo *top_die_info = &die_infos.at(0);
  DieInfo *bottom_die_info = &die_infos.at(1);
  assert(top_die_info->tech_info->lib_cell_num == bottom_die_info->tech_info->lib_cell_num);
  for (int i = 0; i < top_die_info->tech_info->lib_cell_num; ++i) {
    string lib_cell_name;
    int width, height;
    LibCellInfo *lib_cell_info_top = &top_die_info->tech_info->lib_cell_infos.at(i);
    LibCellInfo *lib_cell_info_bottom = &bottom_die_info->tech_info->lib_cell_infos.at(i);
    assert(lib_cell_info_top->name == lib_cell_info_bottom->name);
    lib_cell_name = lib_cell_info_top->name;
    width = floor((lib_cell_info_top->width + lib_cell_info_bottom->width) / 2);
    height = floor((lib_cell_info_top->height + lib_cell_info_bottom->height) / 2);

    dbMaster *master = dbMaster::create(db_lib, lib_cell_name.c_str());
    master->setWidth(width);
    master->setHeight(height);
    master->setType(dbMasterType::CORE);
    // read pins in one LIbCell
    assert(lib_cell_info_top->pin_number == lib_cell_info_bottom->pin_number);
    for (int j = 0; j < lib_cell_info_top->pin_number; ++j) {
      LibPinInfo *pin_info_top = &lib_cell_info_top->lib_pin_infos.at(j);
      LibPinInfo *pin_info_bottom = &lib_cell_info_bottom->lib_pin_infos.at(j);
      assert(pin_info_top->pin_name == pin_info_bottom->pin_name);
      string pin_name = pin_info_top->pin_name;
      int pin_location_x = floor((pin_info_top->pin_location_x + pin_info_bottom->pin_location_x) / 2);
      int pin_location_y = floor((pin_info_top->pin_location_y + pin_info_bottom->pin_location_y) / 2);
      assert(width > pin_location_x);
      assert(height > pin_location_y);

      // (refer to `void lefin::pin` function in submodule/OpenDB/src/lefin/lefin.cpp)
      dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
      dbSigType sig_type = dbSigType::SIGNAL;  // There's no information in this contest benchmarks.
      dbMTerm *master_terminal = dbMTerm::create(master, pin_name.c_str(), io_type, sig_type);
      dbMPin *db_m_pin = dbMPin::create(master_terminal);
      // (refer to `bool lefin::addGeoms` function in submodule/OpenDB/src/lefin/lefin.cpp in case of `lefiGeomRectE`)
      dbBox::create(db_m_pin, db_tech_layer,
                    pin_location_x,
                    pin_location_y,
                    pin_location_x + 1,
                    pin_location_y + 1);
    }
    master->setFrozen();
  }

  // Die Size Setting //
  // for top and bottom die
  pair<int, int> lower_left_pseudo = {0, 0};
  pair<int, int> upper_right_pseudo = {0, 0};

  for (int i = 0; i < 2; ++i) {
    DieInfo &die_info = die_infos.at(i);
    odb::Point lower_left = odb::Point(die_info.lower_left_x, die_info.lower_left_y);
    odb::Point upper_right = odb::Point(die_info.upper_right_x, die_info.upper_right_y);
    odb::Rect die_rect = odb::Rect(lower_left, upper_right);
    db_databases_.at(i)->getChip()->getBlock()->setDieArea(die_rect);

    lower_left_pseudo.first += die_info.lower_left_x;
    lower_left_pseudo.second += die_info.lower_left_y;
    upper_right_pseudo.first += die_info.upper_right_x;
    upper_right_pseudo.second += die_info.upper_right_y;
  }
  lower_left_pseudo.first /= 2;
  lower_left_pseudo.second /= 2;
  upper_right_pseudo.first /= 2;
  upper_right_pseudo.second /= 2;

  // for pseudo die
  odb::Point lower_left_point_pseudo = odb::Point(lower_left_pseudo.first, lower_left_pseudo.second);
  odb::Point upper_right_point_pseudo = odb::Point(upper_right_pseudo.first, upper_right_pseudo.second);
  odb::Rect pseudo_die_rect = odb::Rect(lower_left_point_pseudo, upper_right_point_pseudo);
  db_database_->getChip()->getBlock()->setDieArea(pseudo_die_rect);


  // Instance Setting //
  // for top and bottom die, that will be implemented only when writing lef and def, the end of process
  // for pseudo die
  for (int i = 0; i < instance_number_; ++i) {
    InstanceInfo *instance_info = &instance_infos.at(i);
    dbMaster *master = db_database_->findMaster(instance_info->lib_cell_name.c_str());
    dbInst::create(db_block, master, instance_info->inst_name.c_str());
  }

  // Net and connections setting //
  for (int i = 0; i < net_number_; ++i) {
    // (refer to `dbDatabase* create2LevetDbNoBTerms()` function in submodule/OpenDB/test/cpp/helper.cpp)
    NetInfo *net_info = &net_infos.at(i);
    dbNet *net = dbNet::create(db_block, net_info->net_name.c_str());

    // read pins in one Net
    for (int j = 0; j < net_info->pin_num; ++j) {
      ConnectedPinInfo *pin_info = &net_info->connected_pins_infos.at(j);
      dbInst *inst = db_block->findInst(pin_info->instance_name.c_str());
      inst->findITerm(pin_info->lib_pin_name.c_str())->connect(net);
      assert(inst->findITerm(pin_info->lib_pin_name.c_str()));
    }
  }

  // Row info setting //
  assert(row_infos_.first.row_width == row_infos_.second.row_width);
  assert(row_infos_.first.start_x == row_infos_.second.start_x);
//  assert(row_infos_.first.row_height * row_infos_.first.repeat_count
//             == row_infos_.second.row_height * row_infos_.second.repeat_count);
  // for top and bottom die
  for (int i = 0; i < row_infos_.first.repeat_count; ++i) {
    dbSite *site = dbSite::create(db_databases_.at(0)->findLib("lib"), ("site" + to_string(i)).c_str());
    site->setHeight(row_infos_.first.row_height);
    dbRow::create(db_databases_.at(0)->getChip()->getBlock(), ("row" + to_string(i)).c_str(), site,
                  0, 0, dbOrientType::MX, dbRowDir::HORIZONTAL, 1, row_infos_.first.row_height);
  }
  for (int i = 0; i < row_infos_.second.repeat_count; ++i) {
    dbSite *site = dbSite::create(db_databases_.at(1)->findLib("lib"), ("site" + to_string(i)).c_str());
    site->setHeight(row_infos_.second.row_height);
    dbRow::create(db_databases_.at(1)->getChip()->getBlock(), ("row" + to_string(i)).c_str(), site,
                  0, 0, dbOrientType::MX, dbRowDir::HORIZONTAL, 1, row_infos_.first.row_height);
  }
  // for pseudo die
  int die_height = row_infos_.first.row_height * row_infos_.second.repeat_count;
  int row_height = floor((row_infos_.first.row_height + row_infos_.second.row_height) / 2);
  int repeat_count = floor(die_height / row_height);
  for (int i = 0; i < repeat_count; ++i) {
    dbSite *site = dbSite::create(db_database_->findLib("pseudoDieLib"), ("site" + to_string(i)).c_str());
    site->setHeight(row_height);
    dbRow::create(db_block, ("row" + to_string(i)).c_str(), site,
                  0, 0, dbOrientType::MX, dbRowDir::HORIZONTAL, 1, row_height);
  }

  // Terminal info setting //
  hybrid_size_x_ = terminal_info.size_x;
  hybrid_size_y_ = terminal_info.size_y;
  hybrid_spacing_ = terminal_info.spacing_size;

  // parse end

  init();
}
void Chip::writeICCAD(const string &output_file_name) {
  // open output file
  ofstream ofs(output_file_name);
  assert(ofs.is_open());

  // get the number of instances in the top and bottom die
  int num_instances_top = 0;
  int num_instances_bottom = 0;
  int num_terminals = 0;
  vector<Instance *> terminals;
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    Instance *instance = instance_pointers_.at(i);
    if (instance->getDieId() == 1) {
      num_instances_top++;
    } else if (instance->getDieId() == 2) {
      num_instances_bottom++;
    } else if (instance->isHybridBond()) {
      num_terminals++;
      terminals.push_back(instance);
    }
  }

  // For top die

  ofs << "TopDiePlacement " << num_instances_top << endl;
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    if (instance_pointers_.at(i)->getDieId() != 1)
      continue;
    Instance *instance = instance_pointers_.at(i);
    ofs << "Inst " << instance->getName() << " "
        << instance->getCoordinate().first << " " << instance->getCoordinate().second << endl;
  }

  ofs << "BottomDiePlacement " << num_instances_top << endl;
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    if (instance_pointers_.at(i)->getDieId() != 2)
      continue;
    Instance *instance = instance_pointers_.at(i);
    ofs << "Inst " << instance->getName() << " "
        << instance->getCoordinate().first << " " << instance->getCoordinate().second << endl;
  }

  ofs << "NumTerminals " << num_terminals << endl;
  for (int i = 0; i < num_terminals; ++i) {
    Instance *terminal = terminals.at(i);
    string net_name = terminal->getHybridBondPin()->getNet()->getName();
    auto position = terminal->getHybridBondPin()->getCoordinate();
    ofs << "Terminal " << net_name << " " << position.first << " " << position.second << endl;
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

  i1->findITerm("a")->connect(n1);
  i2->findITerm("b")->connect(n2);
  i3->findITerm("a")->connect(n3);
  i2->findITerm("b")->connect(n4);
  i3->findITerm("a")->connect(n5);
  i3->findITerm("b")->connect(n6);
  i1->findITerm("o")->connect(n5);
  i2->findITerm("o")->connect(n6);
  i3->findITerm("o")->connect(n7);

  i1->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i1->setLocation(30, 20);
  i2->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i2->setLocation(200, 300);
  i3->setPlacementStatus(odb::dbPlacementStatus::PLACED);
  i3->setLocation(0, 0);

  cout << "chip: " << chip << "\t" << db_database_->getChip() << endl;
  cout << "block: " << block << "\t" << db_database_->getChip()->getBlock() << endl;

}
odb::defout::Version Parser::stringToDefVersion(const string &version) {
  if (version == "5.8")
    return odb::defout::Version::DEF_5_8;
  else if (version == "5.7")
    return odb::defout::Version::DEF_5_6;
  else if (version == "5.6")
    return odb::defout::Version::DEF_5_6;
  else if (version == "5.5")
    return odb::defout::Version::DEF_5_5;
  else if (version == "5.4")
    return odb::defout::Version::DEF_5_4;
  else if (version == "5.3")
    return odb::defout::Version::DEF_5_3;
  else
    return odb::defout::Version::DEF_5_8;
}
void Parser::readLef(const string &filename) const {
  odb::lefin lef_reader(db_database_, logger_, false);
  odb::dbLib *lib = lef_reader.createTechAndLib("nangate", filename.c_str());
  odb::dbTech *tech = db_database_->getTech();

  // both are null on parser_ failure
  if (lib != nullptr || tech != nullptr) {
    std::cout << "Lef parsing is succeed." << std::endl;
  } else {
    std::cout << "Lef parsing is failed." << std::endl;
  }
}
void Parser::readDef(const string &filename) const {
  odb::defin def_reader(db_database_, logger_);
  std::vector<odb::dbLib *> search_libs;
  for (odb::dbLib *lib : db_database_->getLibs())
    search_libs.push_back(lib);
  odb::dbChip *chip = def_reader.createChip(search_libs, filename.c_str());
  if (chip) {
    odb::dbBlock *block = chip->getBlock();
    std::cout << "Def parsing is succeed." << std::endl;
  } else {
    std::cout << "Def parsing is failed." << std::endl;
  }
}
void Parser::writeDef(const string &filename, const string &version) const {
  odb::dbChip *chip = db_database_->getChip();
  if (chip) {
    odb::dbBlock *block = chip->getBlock();
    if (block) {
      odb::defout def_writer(logger_);
      def_writer.setVersion(stringToDefVersion(version));
      def_writer.writeBlock(block, filename.c_str());
    }
  } else {
    std::cout << "Writing Def is failed." << std::endl;
  }
}
utl::Logger *Parser::getLogger() const {
  return logger_;
}
void Parser::setLoggerPtr(utl::Logger *logger) {
  logger_ = logger;
}
} // VLSI_backend