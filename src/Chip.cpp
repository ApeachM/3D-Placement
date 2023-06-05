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
#include <algorithm>
#include <fstream>
#include <utility>
#include "Chip.h"
#include "NesterovPlacer.h"
#include "InitialPlacer.h"
#include "Partitioner.h"
#include "HierRTLMP.h"

using namespace std;

namespace VLSI_backend {
Chip::Chip() {
  parser_.setLoggerPtr(&logger_);
}
// main flow //
void Chip::do3DPlace() {

  // 0. target density setting
  // manually setting in code level
  vector<double> densities;
  densities.push_back(2.0); // pseudo die util = top die util + bottom die util
  densities.push_back(1.0);
  densities.push_back(1.0);
  setTargetDensity(densities);

  // 1. do3DPlace the cells in the pseudo die
  phase_ = PHASE::INITIAL_PLACE;
  this->normalPlacement();

/*
  // 2. partition
  phase_ = PHASE::PARTITION;
  this->partition();  // temporary code

  // 3. hybrid bond generate and placement
  phase_ = PHASE::GENERATE_HYBRID_BOND;
  this->generateHybridBonds();
*/

/*
  // 4. placement synchronously
  phase_ = PHASE::TWO_DIE_PLACE;
  this->placement2DieSynchronously();
*/

  phase_ = PHASE::END;
}
void Chip::normalPlacement() {
//  if (!checkDbFile()) {
  doInitialPlace();
  doNesterovPlace();
//  } else {
//    this->dbCaptureRead("db_" + design_name_);
//  }
  this->drawDies();
}
void Chip::partition() {
/*
  partitionSimple();
*/
  if (!checkPartitionFile()) {
    partitioner_ = new Partitioner(nullptr, db_database_, nullptr, &logger_);
    partitioner_->init(design_name_);
    partitioner_->doPartitioning();
    partitioner_->writeSolution();
    delete partitioner_;
  } else {
    readPartitionFile();
  }
/*
  auto *sta = new sta::dbSta;
  sta::dbNetwork *network = sta->getDbNetwork();

  par::PartitionMgr *tritonpart;
  hier_rtl_ = new HierRTLMPartition(network, db_database_, sta, &logger_, tritonpart);
  hier_rtl_->init();
  hier_rtl_->partition();

  delete hier_rtl_;
*/
}
void Chip::generateHybridBonds() {
  // reserve hybrid_bonds and hybrid_bond_pins for preventing to change the addresses
  data_storage_.hybrid_bonds.reserve(net_pointers_.size());

  int hybrid_num = 0;

  // check whether the partitioning was done or not
  for (Instance *instance : instance_pointers_) {
    assert(instance->getDieId() != 0);
  }

  // detect any intersection on the nets
  // (if two dies share one net, then we pretend there is intersection in that net),
  // and if it exists, then make the hybrid bond in the center of the net.
  for (Net *net : net_pointers_) {
    bool intersection = false;
    int die_id_curser = 0;
    for (Pin *pin : net->getConnectedPins()) {
      if (!pin->isInstancePin())
        continue;
      Instance *instance = pin->getInstance();
      net->setDieId(instance->getDieId());
      if (die_id_curser == 0)
        die_id_curser = instance->getDieId();
      else if (die_id_curser != instance->getDieId())
        intersection = true;
    }

    // If there is an intersection in this net,
    // we make a hybrid bond for this net.
    if (intersection) {
      hybrid_num += 1;
      net->setAsIntersected();

      // make objects for hybrid bond
      HybridBond hybrid_bond_object(hybrid_size_x_, hybrid_size_y_, hybrid_spacing_);
      // store them in storage
      data_storage_.hybrid_bonds.push_back(hybrid_bond_object);
      HybridBond *hybrid_bond = &data_storage_.hybrid_bonds.at(data_storage_.hybrid_bonds.size() - 1);

      // set name
      hybrid_bond->setName("hybrid_bond_" + to_string(hybrid_num));

      // link hybrid bond and net
      net->setHybridBond(hybrid_bond); // pin -> net
      hybrid_bond->setConnectedNet(net); // instance -> net

      // set coordinate of hybrid bond
      // p.s. net box would be updated in the first placement phase (Nesterov in virtual die)
      int center_x = static_cast<int>((net->ux() + net->lx()) / 2);
      int center_y = static_cast<int>((net->uy() + net->ly()) / 2);
      hybrid_bond->setCoordinate({center_x, center_y});

      // for iteration
      hybrid_bond_pointers_.push_back(hybrid_bond);
    }
  }

  assert(hybrid_num == data_storage_.hybrid_bonds.size());
  cout << "Hybrid bond #: " << data_storage_.hybrid_bonds.size() << endl;
  this->drawDies("after_hybrid_bond_generation", false, true);
}
void Chip::placement2DieSynchronously() {
  // Now, the initial placement is done by `normalPlacement()`,
  // and the hybrid bond is generated by `generateHybridBonds()`.

  // separate the object pointers for each die
  struct ObjectPointers {
    vector<Instance *> instance_pointers;
    std::vector<Net *> net_pointers_;
    std::vector<Pin *> pin_pointers_;  // This vector includes instance pin pointers and pad pin pointers
    std::vector<Pin *> pad_pointers_;
  };

  ObjectPointers dieVar1, dieVar2;

  for (Instance *instance : instance_pointers_) {
    int die_id = instance->getDieId();
    if (die_id == 1) {
      dieVar1.instance_pointers.push_back(instance);
    } else if (die_id == 2) {
      dieVar2.instance_pointers.push_back(instance);
    }
  }

  for (Net *net : net_pointers_) {
    int die_id = 0;
    if (net->isIntersected()) {
      die_id = -1;
    } else {
      for (Pin *pin : net->getConnectedPins())
        if (pin->isInstancePin()) {
          die_id = pin->getInstance()->getDieId();
          break;
        }
    }
    if (die_id == DIE_ID::INTERSECTED) {
      // the net is on top die and bottom die both
      dieVar1.net_pointers_.push_back(net);
      dieVar2.net_pointers_.push_back(net);
    } else if (die_id == DIE_ID::TOP_DIE)
      // the net is on die top die
      dieVar1.net_pointers_.push_back(net);
    else if (die_id == DIE_ID::BOTTOM_DIE)
      // the net is on bottom die
      dieVar2.net_pointers_.push_back(net);
    else
      // there will be no net on the virtual die on this step.
      assert(0);
  }

  for (Pin *pin : pin_pointers_) {
    int die_id = 0;
    if (!pin->getNet())
      // pass floating pins
      continue;
    die_id = pin->getNet()->getDieId();
    if (die_id == DIE_ID::TOP_DIE) {
      dieVar1.pin_pointers_.push_back(pin);
    } else if (die_id == DIE_ID::BOTTOM_DIE) {
      dieVar2.pin_pointers_.push_back(pin);
    } else if (die_id == DIE_ID::INTERSECTED) {
      dieVar1.pin_pointers_.push_back(pin);
      dieVar2.pin_pointers_.push_back(pin);
    } else
      assert(0);
  }

  // skip the pad pointers because the nesterov optimizer doesn't use them.


  /////////////////////////////////////////////////////////////////////////
  NesterovPlacer nesterov_placer1(
      this->db_database_,
      dieVar1.instance_pointers,
      dieVar1.net_pointers_,
      dieVar1.pin_pointers_,
      dieVar1.pad_pointers_,
      this->die_pointers_.at(DIE_ID::TOP_DIE)
  );
  NesterovPlacer nesterov_placer2(
      this->db_database_,
      dieVar2.instance_pointers,
      dieVar2.net_pointers_,
      dieVar2.pin_pointers_,
      dieVar2.pad_pointers_,
      this->die_pointers_.at(DIE_ID::BOTTOM_DIE)
  );

  nesterov_placer1.setDebugMode(true);
  nesterov_placer2.setDebugMode(true);
  nesterov_placer1.initNesterovPlace(false);
  nesterov_placer2.initNesterovPlace(false);
  nesterov_placer1.updateDB();
  nesterov_placer2.updateDB();
  updateHybridBondPositions();

  if (nesterov_placer1.getMaxNesterovIter() != nesterov_placer2.getMaxNesterovIter())
    assert(0);

  for (int i = 0; i < nesterov_placer1.getMaxNesterovIter(); ++i) {
    int nesterov_iter1, nesterov_iter2;
    nesterov_iter1 = nesterov_placer1.doNesterovPlace(i, true);
    nesterov_iter2 = nesterov_placer2.doNesterovPlace(i, true);
    if (nesterov_iter1 >= nesterov_placer1.getMaxNesterovIter()
        || nesterov_iter2 >= nesterov_placer2.getMaxNesterovIter())
      break;
    nesterov_placer1.updateDB();
    nesterov_placer2.updateDB();
    updateHybridBondPositions();
    cout << "[HPWL]: " << this->getHPWL() << endl;

    string file_name;
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i;
    ss >> file_name;;
    // this->drawDies(file_name, false, true);
  }
  nesterov_placer1.updateDB();
  nesterov_placer2.updateDB();
}

void Chip::init() {
  dbBlock *block = db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();
  dbSet<dbNet> db_nets = block->getNets();

  /**
   * @brief
   * Instance setting
   *
   * @details
   * 1. It makes real instance data and store in \c data_storage.instances. \n
   * 2. Then it makes pointer set for \c Chip class, \n
   * 3. And also makes mapping_ from \c db_instance to instance pointer. \n\n
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
    mapping_.inst_map[instance->getDbInst()] = instance;
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
      mapping_.pin_map_i[pin_pointer->getDbITerm()] = pin_pointer;
    } else if (pin_pointer->isBlockPin()) {
      mapping_.pin_map_b[pin_pointer->getDbBTerm()] = pin_pointer;
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
    mapping_.net_map[net_pointer->getDbNet()] = net_pointer;
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
      Net *net = mapping_.net_map[db_i_term->getNet()];
      Pin *pin = mapping_.pin_map_i[db_i_term];
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
      connected_net = mapping_.net_map[pin->getDbBTerm()->getNet()];
    } else if (pin->isInstancePin()) {
      connected_instance = mapping_.inst_map[pin->getDbITerm()->getInst()];
      connected_net = mapping_.net_map[pin->getDbITerm()->getNet()];
    }
    pin->setConnectedInstance(connected_instance);
    pin->setConnectedNet(connected_net);
  }

  // net -> instances  &&  net -> pins
  for (Net *net : net_pointers_) {
    vector<Instance *> instances;
    vector<Pin *> pins;
    for (dbITerm *i_term : net->getDbNet()->getITerms()) {
      Instance *instance = mapping_.inst_map[i_term->getInst()];
      Pin *pin = mapping_.pin_map_i[i_term];
      instances.push_back(instance);
      pins.push_back(pin);
    }
    for (dbBTerm *b_term : net->getDbNet()->getBTerms()) {
      Pin *pin = mapping_.pin_map_b[b_term];
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

  printDataInfo();
}
ulong Chip::getHPWL() {
  ulong HPWL = 0;
  // TODO: update including intersected die
  for (Net *net : net_pointers_) {
    net->updateBox();
    HPWL += net->getHPWL();
  }
  return HPWL;
}
int Chip::getUnitOfMicro() const {
  return db_database_->getTech()->getDbUnitsPerMicron();
}
void Chip::drawDies(const string &die_name, bool as_dot, bool draw_same_canvas) {
  int scale_factor;
  int die_height_fix = 500;
  if (phase_ < PHASE::GENERATE_HYBRID_BOND) {
    // pseudo die drawing mode
    // let the pixel of the die height be 500
    scale_factor = static_cast<int>(die_pointers_.at(DIE_ID::PSEUDO_DIE)->getHeight() / die_height_fix);
  } else if (phase_ >= PHASE::GENERATE_HYBRID_BOND) {
    // top and bottom die drawing mode
    // let the pixel of the die height be 2000
    // note: top and bottom die size is same here
    scale_factor = static_cast<int>(die_pointers_.at(DIE_ID::TOP_DIE)->getHeight() / die_height_fix);
  }

  if (scale_factor == 0) scale_factor = 10;

  if (phase_ >= PHASE::GENERATE_HYBRID_BOND) {
    uint top_die_w = die_pointers_.at(DIE_ID::TOP_DIE)->getWidth() / scale_factor;
    uint top_die_h = die_pointers_.at(DIE_ID::TOP_DIE)->getHeight() / scale_factor;
    uint bottom_die_w = die_pointers_.at(DIE_ID::BOTTOM_DIE)->getWidth() / scale_factor;
    uint bottom_die_h = die_pointers_.at(DIE_ID::BOTTOM_DIE)->getHeight() / scale_factor;

    Drawer top_die(top_die_w, top_die_h);
    Drawer bottom_die(bottom_die_w, bottom_die_h);

    if (draw_same_canvas) {
      bottom_die.linkImg(top_die.getImage());
      top_die.setCellColor(Color::BLACK);
      bottom_die.setCellColor(Color::RED);
    } else {
      top_die.setCellColor(Color::BLACK);
      bottom_die.setCellColor(Color::BLACK);
    }

    // ID setting
    top_die.setDieId(DIE_ID::TOP_DIE);
    bottom_die.setDieId(DIE_ID::BOTTOM_DIE);

    // Draw cells
    for (Instance *instance : instance_pointers_) {
      int ll_x = instance->getCoordinate().first / scale_factor;
      int ll_y = instance->getCoordinate().second / scale_factor;
      int ur_x = ll_x + instance->getWidth() / scale_factor;
      int ur_y = ll_y + instance->getHeight() / scale_factor;

      if (instance->getDieId() == 1) {
        // top die
        if (as_dot)
          top_die.drawCell(ll_x, ll_y, ll_x, ll_y);
        else
          top_die.drawCell(ll_x, ll_y, ur_x, ur_y);
      } else if (instance->getDieId() == 2) {
        // bottom die
        if (as_dot)
          bottom_die.drawCell(ll_x, ll_y, ll_x, ll_y);
        else
          bottom_die.drawCell(ll_x, ll_y, ur_x, ur_y);
      } else {
        assert(0);
      }
    }
    // Draw hbrid bonds
    for (HybridBond *hybrid_bond : hybrid_bond_pointers_) {
      int ll_x = hybrid_bond->getCoordinate().first / scale_factor;
      int ll_y = hybrid_bond->getCoordinate().second / scale_factor;
      int ur_x = (hybrid_bond->getCoordinate().first + hybrid_size_x_) / scale_factor;
      int ur_y = (hybrid_bond->getCoordinate().second + hybrid_size_y_) / scale_factor;
      if (as_dot)
        top_die.drawHybridBond(ll_x, ll_y, ll_x, ll_y);
      else
        top_die.drawHybridBond(ll_x, ll_y, ur_x, ur_y);
    }

    if (draw_same_canvas)
      top_die.saveImg(die_name);
    else {
      top_die.saveImg("top" + die_name);
      bottom_die.saveImg("bottom" + die_name);
    }
  } else {
    uint pseudo_die_w = die_pointers_.at(DIE_ID::PSEUDO_DIE)->getWidth() / scale_factor;
    uint pseudo_die_h = die_pointers_.at(DIE_ID::PSEUDO_DIE)->getHeight() / scale_factor;

    Drawer pseudo_die(pseudo_die_w, pseudo_die_h);
    pseudo_die.setDieId(0);

    // Draw cells
    for (Instance *instance : instance_pointers_) {
      int ll_x = instance->getCoordinate().first / scale_factor;
      int ll_y = instance->getCoordinate().second / scale_factor;
      int ur_x = ll_x + instance->getWidth() / scale_factor;
      int ur_y = ll_y + instance->getHeight() / scale_factor;
      if (instance->getDieId() != 0)
        assert(0);
      if (as_dot)
        pseudo_die.drawCell(ll_x, ll_y, ll_x + 1, ll_y + 1);
      else
        pseudo_die.drawCell(ll_x, ll_y, ur_x, ur_y);

    }

    pseudo_die.saveImg(die_name);
  }

}
void Chip::printDataInfo() const {
  cout << "======================" << endl;
  cout << "Instance #: " << instance_pointers_.size() << endl;
  cout << "Net #: " << net_pointers_.size() << endl;
  cout << "Pin #:" << pin_pointers_.size() << endl;
  cout << "Die size (x, y): "
       << die_pointers_.at(DIE_ID::PSEUDO_DIE)->getWidth() << ", "
       << die_pointers_.at(DIE_ID::PSEUDO_DIE)->getHeight() << endl;
  cout << "======================" << endl;
}
void Chip::partitionSimple() {
  /* Temporal code */
  int cell_num = static_cast<int>(instance_pointers_.size());
  for (int i = 0; i < floor(cell_num / 2); ++i) {
    Instance *instance = instance_pointers_.at(i);
    instance->assignDie(1);
  }
  for (int i = floor(cell_num / 2); i < cell_num; ++i) {
    Instance *instance = instance_pointers_.at(i);
    instance->assignDie(2);
  }
}
bool Chip::checkPartitionFile() {
  ifstream partition_info_file(design_name_ + "_partition_info");
  if (partition_info_file.fail())
    return false;
  else
    return true;
}
void Chip::readPartitionFile() {
  ifstream partition_info_file(design_name_ + "_partition_info");
  if (partition_info_file.fail())
    assert(0);

  string instance_name;
  int partition_info;
  for (int i = 0; i < instance_number_; ++i) {
    partition_info_file >> instance_name >> partition_info;
    dbInst *db_inst = db_database_->getChip()->getBlock()->findInst(instance_name.c_str());
    Instance *instance = mapping_.inst_map[db_inst];
    instance->assignDie(partition_info + 1);
  }

  // TODO: do for BlockTerminals. This will be for BENCH_TYPE::NORMAL case.
}

// initial placer //
Chip::InitialPlacer::InitialPlacer(std::vector<Instance *> instance_pointers, std::vector<Net *> net_pointers,
                                   std::vector<Pin *> pin_pointers, std::vector<Pin *> pad_pointers,
                                   std::vector<Die *> die_pointers) {
  instance_pointers_ = instance_pointers;
  net_pointers_ = net_pointers;
  pin_pointers_ = pin_pointers;
  pad_pointers_ = pad_pointers;
  die_pointers_ = die_pointers;
}
void Chip::InitialPlacer::placeInstancesCenter() {
  const int center_x = floor(die_pointers_.at(0)->getWidth() / 2);
  const int center_y = floor(die_pointers_.at(0)->getHeight() / 2);

  for (Instance *instance : instance_pointers_) {
    if (!instance->isLocked())
      instance->setCoordinate(center_x - (instance->getWidth() / 2), center_y - (instance->getHeight() / 2));
  }
}
void Chip::InitialPlacer::updatePinInfo() {
  // set the pin is at the boundary of the bounded box or not

  // reset all MinMax attributes
  for (Pin *pin : pin_pointers_) {
    pin->setMinPinXField(false);
    pin->setMinPinYField(false);
    pin->setMaxPinXField(false);
    pin->setMaxPinYField(false);
  }

  for (Net *net : net_pointers_) {
    Pin *pin_min_x = nullptr, *pin_min_y = nullptr;
    Pin *pin_max_x = nullptr, *pin_max_y = nullptr;

    int lx = INT_MAX, ly = INT_MAX;
    int ux = INT_MIN, uy = INT_MIN;

    // Mark B2B info on Pin structures
    for (Pin *pin : net->getConnectedPins()) {
      if (lx > pin->getCoordinate().first) {
        if (pin_min_x)
          pin_min_x->setMinPinXField(false);
        lx = pin->getCoordinate().first;
        pin_min_x = pin;
        pin_min_x->setMinPinXField(true);
      }

      if (ux < pin->getCoordinate().first) {
        if (pin_max_x)
          pin_max_x->setMaxPinXField(false);
        ux = pin->getCoordinate().first;
        pin_max_x = pin;
        pin_max_x->setMaxPinXField(true);
      }

      if (ly > pin->getCoordinate().second) {
        if (pin_min_y)
          pin_min_y->setMinPinYField(false);
        ly = pin->getCoordinate().second;
        pin_min_y = pin;
        pin_min_y->setMinPinYField(true);
      }

      if (uy < pin->getCoordinate().second) {
        if (pin_max_y)
          pin_max_y->setMaxPinYField(false);
        uy = pin->getCoordinate().second;
        pin_max_y = pin;
        pin_max_y->setMaxPinYField(true);
      }
    }
  }
}
void Chip::InitialPlacer::createSparseMatrix() {
  // This function is from below link
  // https://github.com/The-OpenROAD-Project/OpenROAD/blob/977c0794af50e0d3ed993d324b0adead87e32782/src/gpl/src/initialPlace.cpp#L234
  const int placeCnt = instance_pointers_.size();
  instLocVecX_.resize(placeCnt);
  fixedInstForceVecX_.resize(placeCnt);
  instLocVecY_.resize(placeCnt);
  fixedInstForceVecY_.resize(placeCnt);

  placeInstForceMatrixX_.resize(placeCnt, placeCnt);
  placeInstForceMatrixY_.resize(placeCnt, placeCnt);

  //
  // list_x and list_y is a temporary vector that have tuples, (idx1, idx2, val)
  //
  // list_x finally becomes placeInstForceMatrixX_
  // list_y finally becomes placeInstForceMatrixY_
  //
  // The triplet vector is recommended usages
  // to fill in SparseMatrix from Eigen docs.
  //
  vector<T> list_x, list_y;
  list_x.reserve(1000000);
  list_y.reserve(1000000);

  // initialize vector
  for (Instance *instance : instance_pointers_) {
    int idx = instance->getId();
    instLocVecX_(idx) = instance->getCenterX();
    instLocVecY_(idx) = instance->getCenterY();

    fixedInstForceVecX_(idx) = 0;
    fixedInstForceVecY_(idx) = 0;
  }

  // for each net
  for (Net *net : net_pointers_) {
    // skip for small nets.
    if (net->getConnectedPins().size() <= 1)
      continue;

    // escape long time cals on huge fanout.
    if (net->getConnectedPins().size() >= max_fan_out_)
      continue;

    float net_weight = net_weight_scale_ / static_cast<float>(net->getConnectedPins().size() - 1);

    vector<Pin *> pins = net->getConnectedPins();
    for (int pin_idx1 = 1; pin_idx1 < pins.size(); ++pin_idx1) {
      Pin *pin1 = pins.at(pin_idx1);
      for (int pin_idx2 = 0; pin_idx2 < pin_idx1; ++pin_idx2) {
        Pin *pin2 = pins.at(pin_idx2);

        // no need to fill in when instance is same
        if (pin1->getInstance() == pin2->getInstance())
          continue;

        // B2B modeling on min_x/max_x pins.
        if (pin1->isMinPinX() || pin1->isMaxPinX() || pin2->isMinPinX() || pin2->isMaxPinX()) {
          int diff_x = abs(pin1->getCoordinate().first - pin2->getCoordinate().first);
          float weight_x = 0;
          if (diff_x > min_diff_length_)
            weight_x = net_weight / diff_x;
          else
            weight_x = net_weight / min_diff_length_;

          // both pin came from instance
          if (pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            const int inst2 = pin2->getInstance()->getId();

            list_x.push_back(T(inst1, inst2, weight_x));
            list_x.push_back(T(inst2, inst1, weight_x));
            list_x.push_back(T(inst1, inst2, -weight_x));
            list_x.push_back(T(inst2, inst1, -weight_x));

            fixedInstForceVecX_(inst1) +=
                -weight_x * static_cast<float>
                (
                    (pin1->getCoordinate().first - pin1->getInstance()->getCenterX())
                        - (pin2->getCoordinate().first - pin2->getInstance()->getCenterX())
                );
            fixedInstForceVecX_(inst2) +=
                -weight_x * static_cast<float>
                (
                    (pin2->getCoordinate().first - pin2->getInstance()->getCenterX())
                        - (pin1->getCoordinate().first - pin1->getInstance()->getCenterX())
                );
          }
            // pin1 from IO port / pin2 from Instance
          else if (!pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst2 = pin2->getInstance()->getId();
            list_x.push_back(T(inst2, inst2, weight_x));
            fixedInstForceVecX_(inst2) += weight_x * static_cast<float>
            (pin1->getCoordinate().first - (pin2->getCoordinate().first - pin2->getInstance()->getCenterX()));
          }
            // pin1 from Instance / pin2 from IO port
          else if (pin1->isInstancePin() && !pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            list_x.push_back(T(inst1, inst1, weight_x));
            fixedInstForceVecX_(inst1) += weight_x * static_cast<float>
            (pin2->getCoordinate().first - (pin1->getCoordinate().first - pin1->getInstance()->getCenterX()));
          }
        }

        // B2B modeling on min_y/max_y pins.
        if (pin1->isMinPinY() || pin1->isMaxPinY() || pin2->isMinPinY() || pin2->isMaxPinY()) {
          int diff_y = abs(pin1->getCoordinate().second - pin2->getCoordinate().second);
          float weight_y = 0;
          if (diff_y > min_diff_length_) {
            weight_y = net_weight / static_cast<float>(diff_y);
          } else {
            weight_y = net_weight / min_diff_length_;
          }

          // both pin came from instance
          if (pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            const int inst2 = pin2->getInstance()->getId();
            list_y.push_back(T(inst1, inst1, weight_y));
            list_y.push_back(T(inst2, inst2, weight_y));
            list_y.push_back(T(inst1, inst2, -weight_y));
            list_y.push_back(T(inst2, inst1, -weight_y));

            fixedInstForceVecY_(inst1) += -weight_y * static_cast<float>
            (
                (pin1->getCoordinate().second - pin1->getInstance()->getCenterY())
                    - (pin2->getCoordinate().second - pin2->getInstance()->getCenterY())
            );
            fixedInstForceVecY_(inst2) += -weight_y * static_cast<float>
            (
                (pin2->getCoordinate().second - pin2->getInstance()->getCenterY())
                    - (pin1->getCoordinate().second - pin1->getInstance()->getCenterY())
            );
          }
            // pin1 from IO port // pin2 from Instance
          else if (!pin1->isInstancePin() && pin2->isInstancePin()) {
            const int inst2 = pin2->getInstance()->getId();
            list_y.push_back(T(inst2, inst2, weight_y));
            fixedInstForceVecY_(inst2) += weight_y * static_cast<float>(
                (pin1->getCoordinate().second - (pin2->getCoordinate().second - pin2->getInstance()->getCenterY()))
            );
          }
            // pin1 from Instance / pin2 from IO port
          else if (pin1->isInstancePin() && !pin2->isInstancePin()) {
            const int inst1 = pin1->getInstance()->getId();
            list_y.push_back(T(inst1, inst1, weight_y));
            fixedInstForceVecY_(inst1) += weight_y * static_cast<float>(
                (pin2->getCoordinate().second - (pin1->getCoordinate().second - pin1->getInstance()->getCenterY()))
            );

          }
        }
      }
    }
  }
  placeInstForceMatrixX_.setFromTriplets(list_x.begin(), list_x.end());
  placeInstForceMatrixY_.setFromTriplets(list_y.begin(), list_y.end());
}
void Chip::InitialPlacer::setPlaceIDs() {
  // reset ExtId for all instances
  for (Instance *instance : instance_pointers_) {
    instance->setId(INT_MAX);
  }
  // set index only with place-able instances
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    Instance *instance = instance_pointers_.at(i);
    if (!instance->isFixed())
      instance->setId(i);
  }
}
pair<float, float> Chip::InitialPlacer::cpuSparseSolve() {
  pair<float, float> error;
  BiCGSTAB<SMatrix, IdentityPreconditioner> solver;
  solver.setMaxIterations(max_solver_iter_);
  // for x
  solver.compute(placeInstForceMatrixX_);
  instLocVecX_ = solver.solveWithGuess(fixedInstForceVecX_, instLocVecX_);
  error.first = solver.error();
  // for y
  solver.compute(placeInstForceMatrixY_);
  instLocVecY_ = solver.solveWithGuess(fixedInstForceVecY_, instLocVecY_);
  error.second = solver.error();

  return error;
}

// parser //
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
  bench_type_ = BENCH_TYPE::ICCAD;
  setDesignName(input_file_name);

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
        // dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
        // for partitioning, we should make any flow of IO.
        // let the last pin has the output pin, and the others has input flow
        dbIoType io_type;
        if (j != lib_cell_info.pin_number - 1)
          io_type = dbIoType::INPUT;
        else
          io_type = dbIoType::OUTPUT;

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
      // dbIoType io_type = dbIoType::INOUT;  // There's no information in this contest benchmarks.
      // for partitioning, we should make any flow of IO.
      // let the last pin has the output pin, and the others has input flow
      dbIoType io_type;
      if (j != lib_cell_info_top->pin_number - 1)
        io_type = dbIoType::INPUT;
      else
        io_type = dbIoType::OUTPUT;
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
  vector<HybridBond *> terminals;
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    Instance *instance = instance_pointers_.at(i);
    if (instance->getDieId() == 1) {
      num_instances_top++;
    } else if (instance->getDieId() == 2) {
      num_instances_bottom++;
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

  ofs << "NumTerminals " << hybrid_bond_pointers_.size() << endl;
  for (int i = 0; i < hybrid_bond_pointers_.size(); ++i) {
    HybridBond *terminal = terminals.at(i);
    string net_name = terminal->getConnectedNet()->getName();
    auto position = terminal->getCoordinate();
    ofs << "Terminal " << net_name << " " << position.first << " " << position.second << endl;
  }

}
void Chip::test() {
  design_name_ = "db_outputCASE2";
  this->dbCaptureRead(design_name_);
  this->drawDies("test", false, true);;
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

// normal placement //
void Chip::setTargetDensity(vector<double> densities) {
  if (densities.size() != die_pointers_.size())
    assert(0);
  for (int i = 0; i < densities.size(); ++i) {
    die_pointers_.at(i)->setDensity(densities.at(i));
  }
}
void Chip::doInitialPlace() {
  // This function is from below link
  // https://github.com/The-OpenROAD-Project/OpenROAD/blob/977c0794af50e0d3ed993d324b0adead87e32782/src/gpl/src/initialPlace.cpp#L82
  InitialPlacer initial_placer(
      this->instance_pointers_,
      this->net_pointers_,
      this->pin_pointers_,
      this->pad_pointers_,
      this->die_pointers_);

  initial_placer.placeInstancesCenter();
  initial_placer.setPlaceIDs();
  for (int iter = 0; iter < initial_placer.max_iter_; ++iter) {
    initial_placer.updatePinInfo();
    initial_placer.createSparseMatrix();
    pair<float, float> error = initial_placer.cpuSparseSolve();
    float error_max = max(error.first, error.second);
    cout << "[InitialPlace] Iter: " << iter << "\tHPWL: " << getHPWL() << " CG residual: " << error_max << endl;
    if (error_max < 1e-5 && iter >= 5)
      break;
  }
}
void Chip::doNesterovPlace() {
  NesterovPlacer nesterov_placer(
      this->db_database_,
      this->instance_pointers_,
      this->net_pointers_,
      this->pin_pointers_,
      this->pad_pointers_,
      this->die_pointers_.at(DIE_ID::PSEUDO_DIE));
  nesterov_placer.setDebugMode(true);
  nesterov_placer.initNesterovPlace();
  nesterov_placer.setMaxNesterovIter(300);
  nesterov_placer.doNesterovPlace();
  cout << "[HPWL] : " << getHPWL() << endl;
}
int Chip::getInstanceNumber() const {
  return instance_number_;
}
void Chip::setInstanceNumber(int instance_number) {
  data_storage_.instances.reserve(instance_number);
  instance_number_ = instance_number;
}
int Chip::getNetNumber() const {
  return net_number_;
}
void Chip::setNetNumber(int net_number) {
  data_storage_.nets.reserve(net_number);
  net_number_ = net_number;
}
dbDatabase *Chip::getDbDatabase() const {
  return db_database_;
}
void Chip::setDbDatabase(dbDatabase *db_database) {
  parser_.db_database_ = db_database;
  db_database_ = db_database;
}
void Chip::setDesignName(const string &input_file_name) {
  if (bench_type_ == BENCH_TYPE::ICCAD) {
    if (input_file_name == "../test/benchmarks/iccad2022/case2.txt")
      design_name_ = "CASE2";
    else if (input_file_name == "../test/benchmarks/iccad2022/case2_hidden.txt")
      design_name_ = "CASE2_HIDDEN";
    else if (input_file_name == "../test/benchmarks/iccad2022/case3.txt")
      design_name_ = "CASE3";
    else if (input_file_name == "../test/benchmarks/iccad2022/case3_hidden.txt")
      design_name_ = "CASE3_HIDDEN";
    else if (input_file_name == "../test/benchmarks/iccad2022/case4.txt")
      design_name_ = "CASE4";
    else if (input_file_name == "../test/benchmarks/iccad2022/case4_hidden.txt")
      design_name_ = "CASE4_HIDDEN";
  } else {
    design_name_ = "NORMAL"; // TODO: re-specify this
  }
}
void Chip::updateHybridBondPositions() {
  for (HybridBond *hybrid_bond : hybrid_bond_pointers_) {
    hybrid_bond->updatePosition();
  }
}

// Nesterov Placer //
Chip::NesterovPlacer::NesterovPlacer(odb::dbDatabase *db_database,
                                     std::vector<Instance *> instance_pointers,
                                     std::vector<Net *> net_pointers,
                                     std::vector<Pin *> pin_pointers,
                                     std::vector<Pin *> pad_pointers,
                                     Die *die_pointer) {
  db_database_ = db_database;
  // TODO: should check move whether method is proper or not.
  instance_pointers_ = std::move(instance_pointers);
  net_pointers_ = std::move(net_pointers);
  pin_pointers_ = std::move(pin_pointers);
  pad_pointers_ = std::move(pad_pointers);
  die_pointer_ = die_pointer;


  // hyper parameters
  if (instance_pointers_.size() < 1e5)
    initDensityPenalty = 0.01;
  else if (instance_pointers.size() < 1e3)
    initDensityPenalty = 0.1;

}
bool Chip::NesterovPlacer::initNesterovPlace(bool is_pseudo_die) {
  if (!is_base_initialized_) {
    // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L1054
    // bool Replace::initNesterovPlace()
    // void NesterovBase::init()
    is_base_initialized_ = true;

    if (instance_pointers_.empty()) {
      cout << "No placeable instance - skipping placement." << endl;
      return false;
    }
    setInstancesArea();

    // void NesterovBase::init()
    // Set a fixed seed
    srand(42);
    int dbu_per_micron = db_database_->getChip()->getBlock()->getDbUnitsPerMicron();

    if (is_pseudo_die) {
      for (Instance *instance : instance_pointers_) {
        // For any cell, add a random noise between -1 and 1 microns to each of its
        // x and y components. This is added to make it very unlikely that identical
        // cells connected in parallel do not start at the exact same position and
        // consequently shadow each other throughout the entire placement process
        int x_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
        int y_offset = rand() % (2 * dbu_per_micron) - dbu_per_micron;
        instance->setCoordinate(instance->getCoordinate().first + x_offset,
                                instance->getCoordinate().second + y_offset);
      }
    }
    initFillerCells();
    initBins();

    // initialize fft structure based on bins
    fft_ = new gpl::FFT(bin_cnt_x_, bin_cnt_y_, bin_size_x_, bin_size_y_);

    for (Instance *instance : instance_pointers_) {
      instance->setDensityValueAsDefault();
    }
    for (Pin *pin : pin_pointers_) {
      pin->initDensityCoordinate();
    }
    updateDensitySize();
  }

  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L138
  // void NesterovPlace::init()

  initSLPStepsVars();

  for (int idx = 0; idx < instance_pointers_.size(); ++idx) {
    Instance *cell = instance_pointers_.at(idx);
    updateDensityCoordiLayoutInside(cell);
    cur_slp_coordinates_[idx] = prev_slp_coordinates_[idx] = cur_coordinates_[idx] = init_coordinates_[idx] =
        pair<float, float>{cell->getDensityCenterX(), cell->getDensityCenterY()};
  }
  // bin
  updateGCellDensityCenterLocation(cur_slp_coordinates_);
  prev_hpwl_ = getHpwl();
  // FFT update
  updateDensityForceBin();
  base_wire_length_coefficient_ = initWireLengthCoef / (static_cast<float>(bin_size_x_ + bin_size_y_) * 0.5);
  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());
  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());
  updateWireLengthCoef(sum_overflow_);
  // TODO: check whether it is proper to print these in this order, also below one.
  cout << "[replace] np init: InitialHPWL: " << prev_hpwl_ << endl;
  cout << "[replace] np init: BaseWireLengthCoef: " << base_wire_length_coefficient_ << endl;
  cout << "[replace] np init: OverflowArea: " << overflow_area_ << endl;
  cout << "[replace] np init: NesterovInstArea: " << nesterovInstsArea() << endl;
  cout << "[replace] np init: InitSumOverflow: " << sum_overflow_unscaled_ << endl;
  cout << "[replace] np init: WireLengthCoef: " << wire_length_coefficient_x_ << endl;

  // WL update
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  // fill in cur_slp_sum_grads_, cur_slp_wire_length_grads_, cur_slp_density_grads_
  updateGradients(cur_slp_sum_grads_, cur_slp_wire_length_grads_, cur_slp_density_grads_);

  if (is_diverged_) { return false; }

  // approximately fill in prev_slp_coordinates_ to calculate lc vars
  updateInitialPrevSLPCoordi();

  // bin, FFT, wlen update with prevSLPCoordi.
  // prev_slp_coordinates_ 얘 안 바뀌였나..?
  updateGCellDensityCenterLocation(prev_slp_coordinates_);
  updateDensityForceBin();
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  // update previSumGrads_, prev_slp_wire_length_grads_, prev_slp_density_grads_
  updateGradients(prev_slp_sum_grads_, prev_slp_wire_length_grads_, prev_slp_density_grads_);

  if (is_diverged_) { return false; }

  density_penalty_ = (wire_length_grad_sum_ / density_grad_sum_) * initDensityPenalty;
  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());
  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());
  step_length_ = getStepLength();
  cout << "[replace] np init: WireLengthGradSum: " << wire_length_grad_sum_ << endl;
  cout << "[replace] np init: DensityGradSum: " << density_grad_sum_ << endl;
  cout << "[replace] np init: InitDensityPenalty: " << density_penalty_ << endl;
  cout << "[replace] np init: PrevSumOverflow: " << sum_overflow_unscaled_ << endl;
  cout << "[replace] np init: InitialStepLength: " << step_length_ << endl;

  if ((isnan(step_length_) || isinf(step_length_)) && recursion_cnt_init_slp_coef_ < maxRecursionInitSLPCoef) {
    initialPrevCoordiUpdateCoef *= 10;
    cout << "np init: steplength = 0 detected. Rerunning Nesterov::init() with initPrevSLPCoef: "
         << initialPrevCoordiUpdateCoef << endl;
    recursion_cnt_init_slp_coef_++;
    initNesterovPlace(is_pseudo_die);
  }

  if (isnan(step_length_) || isinf(step_length_)) {
    cout << "RePlAce diverged at initial iteration. Re-run with a smaller init_density_penalty value." << endl;
  }

  return true;
}
int Chip::NesterovPlacer::doNesterovPlace(int start_iter, bool only_one_iter) {
  // refer: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/nesterovPlace.cpp#L482
  // int NesterovPlace::doNesterovPlace(int start_iter)

  // backTracking variable.
  // divergence detection
  float minSumOverflow = 1e30;
  float hpwlWithMinSumOverflow = 1e30;

  // dynamic adjustment of max_phi_coef
  bool is_max_phi_coef_changed = false;

  // snapshot saving detection
  bool isSnapshotSaved = false;

  // snapshot info
  vector<pair<float, float>> snapshot_coordinates;
  vector<pair<float, float>> snapshot_slp_coordinates;
  vector<pair<float, float>> snapshot_slp_sum_grads;

  float snapshot_a = 0;
  float snapshot_density_penalty = 0;
  float snapshot_step_length = 0;
  float snapshot_wl_coef_x = 0, snapshot_wl_coef_y = 0;

  bool is_diverge_tried_revert = false;

  // Core Nesterov Loop
  int iter = start_iter;
  for (; iter < max_nesterov_iter_; ++iter) {
    // cout << "[replace-test] np: InitSumOverflow: " << sum_overflow_unscaled_ << endl;

    prevA = curA;
    // here, prevA is a_(k), curA is a_(k+1)
    // See, the ePlace-MS paper's Algorithm 1
    curA = (1.0 + sqrt(4.0 * prevA * prevA + 1.0)) * 0.5;
    // coeff is (a_k - 1) / ( a_(k+1) ) in paper.
    float coeff = (prevA - 1.0) / curA;

    // Back-Tracking loop
    int num_back_trak;
    for (num_back_trak = 0; num_back_trak < max_back_track_; ++num_back_trak) {
      // fill in nextCoordinates with given step_length_
      // here, the instance_pointers_ includes the fillers
      for (int k = 0; k < instance_pointers_.size(); ++k) {
        pair<float, float> next_coordinate(
            cur_slp_coordinates_[k].first + step_length_ * cur_slp_sum_grads_[k].first,
            cur_slp_coordinates_[k].second + step_length_ * cur_slp_sum_grads_[k].second);
        pair<float, float> next_slp_coordinate(
            next_coordinate.first + coeff * (next_coordinate.first - cur_coordinates_[k].first),
            next_coordinate.second + coeff * (next_coordinate.second - cur_coordinates_[k].second));
        Instance *current_instance = instance_pointers_.at(k);

        next_coordinates_.at(k) = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, next_coordinate.first),
            getDensityCoordiLayoutInsideY(current_instance, next_coordinate.second));
        next_slp_coordinates_[k] = pair<float, float>(
            getDensityCoordiLayoutInsideX(current_instance, next_slp_coordinate.first),
            getDensityCoordiLayoutInsideY(current_instance, next_slp_coordinate.second));
      }
      updateGCellDensityCenterLocation(next_slp_coordinates_);
      updateDensityForceBin();
      updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);
      updateGradients(next_slp_sum_grads_, next_slp_wire_length_grads_, next_slp_density_grads_);
      if (is_diverged_)
        break;
      if (stepLengthDivergeCheck())
        break;
    }

    // dynamic adjustment for
    // better convergence with
    // large designs
    if (!is_max_phi_coef_changed && sum_overflow_unscaled_ < 0.35f) {
      is_max_phi_coef_changed = true;
      max_phi_coef_ *= 0.99;
    }
    if (max_back_track_ == num_back_trak) {
      cout << "Backtracking limit reached so a small step will be taken" << endl;
    }
    if (is_diverged_) { break; }

    updateNextIter();
    printStateNesterov(iter);

    if (minSumOverflow > sum_overflow_unscaled_) {
      minSumOverflow = sum_overflow_unscaled_;
      hpwlWithMinSumOverflow = prev_hpwl_;
    }
    /*
     diverge detection on
     large max_phi_cof value + large design

     1) happen overflow < 20%
     2) Hpwl is growing
    */
    if (sum_overflow_unscaled_ < 0.3f && sum_overflow_unscaled_ - minSumOverflow >= 0.02f
        && hpwlWithMinSumOverflow * 1.2f < prev_hpwl_) {
      handleDiverge(snapshot_coordinates, snapshot_slp_coordinates, snapshot_slp_sum_grads, snapshot_a,
                    snapshot_density_penalty, snapshot_step_length, snapshot_wl_coef_x, snapshot_wl_coef_y,
                    is_diverge_tried_revert);
    }
    // if it reached target overflow
    if (finishCheck()) {
      iter = max_nesterov_iter_;
      break;
    }
    if (debug_mode_) {
      string file_name;
      std::stringstream ss;
      ss << std::setw(4) << std::setfill('0') << iter;
      ss >> file_name;;
      if (this->die_pointer_->getDieId() == DIE_ID::TOP_DIE)
        file_name = "top_" + file_name;
      else if (this->die_pointer_->getDieId() == DIE_ID::BOTTOM_DIE)
        file_name = "bottom_" + file_name;
      else if (this->die_pointer_->getDieId() == DIE_ID::PSEUDO_DIE)
        file_name = "pseudo_" + file_name;
      drawDie(file_name);
    }

    if (only_one_iter)
      return iter;
  }
  // in all case including diverge,
  // db should be updated.
  updateDB();
  if (is_diverged_) {
    cout << "log_->error(GPL, diverge_code_, diverge_msg_);" << endl;
  }
  return iter;

}
bool Chip::NesterovPlacer::finishCheck() const {
  if (sum_overflow_unscaled_ <= targetOverflow) {
    cout << "[NesterovSolve] Finished with Overflow: " << sum_overflow_unscaled_ << endl;
    return true;
  }
  return false;
}
void Chip::NesterovPlacer::printStateNesterov(int iter) const {
  if (iter == 0 || (iter + 1) % 10 == 0) {
    cout << "[NesterovSolve] Iter: " << iter + 1
         << "\toverflow: " << sum_overflow_unscaled_
         << "\tHPWL: " << prev_hpwl_
         << endl;
  }
}
bool Chip::NesterovPlacer::stepLengthDivergeCheck() {
  float newStepLength = getStepLength();

  if (isnan(newStepLength) || isinf(newStepLength)) {
    is_diverged_ = true;
    diverge_msg_ = "RePlAce diverged at newStepLength.";
    diverge_code_ = 305;
    return true;
  }
  if (newStepLength > step_length_ * 0.95) {
    step_length_ = newStepLength;
    return true;
  } else if (newStepLength < 0.01) {
    step_length_ = 0.01;
    return true;
  } else {
    step_length_ = newStepLength;
    return false;
  }
}
void Chip::NesterovPlacer::handleDiverge(const vector<pair<float, float>> &snapshotCoordi,
                                         const vector<pair<float, float>> &snapshotSLPCoordi,
                                         const vector<pair<float, float>> &snapshotSLPSumGrads,
                                         float snapshotA,
                                         float snapshotDensityPenalty,
                                         float snapshotStepLength,
                                         float snapshotWlCoefX,
                                         float snapshotWlCoefY,
                                         bool &isDivergeTriedRevert) {
  diverge_msg_ = "RePlAce divergence detected. ";
  diverge_msg_ += "Re-run with a smaller max_phi_cof value.";
  diverge_code_ = 307;
  is_diverged_ = true;


  // revert back the current density penality
  cur_coordinates_ = snapshotCoordi;
  cur_slp_coordinates_ = snapshotSLPCoordi;
  cur_slp_sum_grads_ = snapshotSLPSumGrads;
  curA = snapshotA;
  density_penalty_ = snapshotDensityPenalty;
  step_length_ = snapshotStepLength;
  wire_length_coefficient_x_ = snapshotWlCoefX;
  wire_length_coefficient_y_ = snapshotWlCoefY;

  updateGCellDensityCenterLocation(cur_coordinates_);
  updateDensityForceBin();
  updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

  is_diverged_ = false;
  diverge_code_ = 0;
  diverge_msg_ = "";
  isDivergeTriedRevert = true;

  // turn off the RD forcely
  is_routability_need_ = false;
}

void Chip::NesterovPlacer::setInstancesArea() {
  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/a5e786eb65f40abfb7004b18312d519dac95cc33/src/gpl/src/placerBase.cpp#L798
  // void PlacerBase::init()
  for (Instance *instance : instance_pointers_) {
    if (instance->isInstance()) {
      if (instance->isFixed()) {
        // TODO
        //  we didn't specify the each die for nesterov.
        //  Therefore, we remain this part as blank
        // Check whether fixed instance is
        // within the core area
        // outside core area is none of RePlAce's business
        /*
        if(isCoreAreaOverlap){
          // ... //
        }
        */
      } else {
        uint inst_area = instance->getArea();
        place_instances_area_ += inst_area;
        // macro cells should be
        // macro_instances_area_
        uint row_height = (*db_database_->getChip()->getBlock()->getRows().begin())->getSite()->getHeight();
        if (instance->getHeight() > (row_height * 6)) {
          macro_instances_area_ += inst_area;
        }
          // smaller or equal height cells should be
          // stdInstArea_
        else {
          std_instances_area_ += inst_area;
        }
      }
    } else if (instance->isDummy()) {
      dummyInsts_.push_back(instance);
      nonPlaceInsts_.push_back(instance);
      non_place_instances_area_ += instance->getArea();
    }
  }
}
void Chip::NesterovPlacer::initFillerCells() {
  // extract average dx/dy in range (10%, 90%)
  vector<int> width_storage;
  vector<int> height_storage;
  width_storage.reserve(instance_pointers_.size());
  height_storage.reserve(instance_pointers_.size());
  for (Instance *instance : instance_pointers_) {
    width_storage.push_back(static_cast<int>(instance->getWidth()));
    height_storage.push_back(static_cast<int>(instance->getHeight()));
  }

  // sort
  std::sort(width_storage.begin(), width_storage.end());
  std::sort(height_storage.begin(), height_storage.end());

  // average from (10 - 90%)
  int64_t width_sum = 0, height_sum = 0;

  int min_idx = static_cast<int>(static_cast<float>(width_storage.size()) * 0.05);
  int max_idx = static_cast<int>(static_cast<float>(width_storage.size()) * 0.95);

  // when #instances are too small,
  // extracts average values in whole ranges.
  if (min_idx == max_idx) {
    min_idx = 0;
    max_idx = static_cast<int>(width_storage.size());
  }

  for (int i = min_idx; i < max_idx; i++) {
    width_sum += width_storage[i];
    height_sum += height_storage[i];
  }

  // the avg width and avg height will be used as filler cells
  // width and height
  filler_width_ = static_cast<int>(width_sum / (max_idx - min_idx));
  filler_height_ = static_cast<int>(height_sum / (max_idx - min_idx));

  int64_t core_area = die_pointer_->getArea();

  // nonPlaceInstsArea should not have density downscaling!!!
  white_space_area_ = core_area - non_place_instances_area_;

  if (use_uniform_target_density_) {
    target_density_ =
        static_cast<float>(std_instances_area_) / static_cast<float>(white_space_area_ - macro_instances_area_)
            + 0.01;
  } else {
    target_density_ = die_pointer_->getDensity();
  }

  // density screening
  movable_area_ = white_space_area_ * target_density_;
  total_filler_area_ = movable_area_ - nesterovInstsArea();
  uniform_target_density_ = static_cast<float>(nesterovInstsArea()) / static_cast<float>(white_space_area_);
  if (total_filler_area_ < 0) {
    uniform_target_density_ = ceilf(uniform_target_density_ * 100) / 100;
    cout << "Use a higher -density or" << endl <<
         "re-floorplan with a larger core area." << endl <<
         "Given target density: " << target_density_ << endl <<
         "Suggested target density: " << uniform_target_density_ << endl;
  }
  int fillerCnt = static_cast<int>(
      total_filler_area_ / static_cast<int64_t>(filler_width_ * filler_height_));

  mt19937 genX;
  mt19937 genY;
  if (seed_fix) {
    genX.seed(1111); // fix seed
    genY.seed(2222); // fix seed
  } else {
    genX.seed(random_device{}());
    genY.seed(random_device{}());
  }
  uniform_int_distribution<int> disX(0, (int) die_pointer_->getWidth());
  uniform_int_distribution<int> disY(0, (int) die_pointer_->getHeight());

  // make and store the fillers
  fillers_.reserve(fillerCnt);
  for (int i = 0; i < fillerCnt; ++i) {
    Instance filler;
    filler.setWidth(filler_width_);
    filler.setHeight(filler_height_);
    filler.setCoordinate(disX(genX), disY(genY));
    fillers_.push_back(filler);
  }
  // setting the pointers of fillers in real data storage
  for (auto &filler : fillers_) {
    instance_pointers_.push_back(&filler);
  }

}
void Chip::NesterovPlacer::initBins() {
  // refer to: https://github.com/The-OpenROAD-Project/OpenROAD/blob/402c5cff5d5dac9868f812fec69edb064a5bfbb3/src/gpl/src/nesterovBase.cpp#L723
  // void BinGrid::initBins()
  int ux_ = die_pointer_->getUpperRightX();
  int lx_ = die_pointer_->getLowerLeftX();
  int uy_ = die_pointer_->getUpperRightY();
  int ly_ = die_pointer_->getLowerLeftY();
  int64_t totalBinArea = static_cast<int64_t>(ux_ - lx_) * static_cast<int64_t>(uy_ - ly_);
  int64_t averagePlaceInstArea = place_instances_area_ / (instance_pointers_.size() - fillers_.size());

  int64_t idealBinArea = std::round(static_cast<float>(averagePlaceInstArea) / target_density_);
  int idealBinCnt = totalBinArea / idealBinArea;
  if (idealBinCnt < 4)   // the smallest we allow is 2x2 bins
    idealBinCnt = 4;
/*
      log_->info(GPL, 23, "TargetDensity: {:.2f}", targetDensity_);
      log_->info(GPL, 24, "AveragePlaceInstArea: {}", averagePlaceInstArea);
      log_->info(GPL, 25, "IdealBinArea: {}", idealBinArea);
      log_->info(GPL, 26, "IdealBinCnt: {}", idealBinCnt);
      log_->info(GPL, 27, "TotalBinArea: {}", totalBinArea);
*/
  int foundBinCnt = 2;
  // find binCnt: 2, 4, 8, 16, 32, 64, ...
  // s.t. binCnt^2 <= idealBinCnt <= (binCnt*2)^2.
  for (foundBinCnt = 2; foundBinCnt <= 1024; foundBinCnt *= 2) {
    if (foundBinCnt * foundBinCnt <= idealBinCnt
        && 4 * foundBinCnt * foundBinCnt > idealBinCnt) {
      break;
    }
  }
  bin_cnt_x_ = foundBinCnt;
  bin_cnt_y_ = foundBinCnt;
  // log_->info(GPL, 28, "BinCnt: {} {}", bin_cnt_x_, bin_cnt_y_);
  bin_size_x_ = ceil(static_cast<float>((ux_ - lx_)) / bin_cnt_x_);
  bin_size_y_ = ceil(static_cast<float>((uy_ - ly_)) / bin_cnt_y_);
  // log_->info(GPL, 29, "BinSize: {} {}", bin_size_x_, bin_size_y_);
  binStor_.resize(bin_cnt_x_ * bin_cnt_y_);
  bins_.reserve(bin_cnt_x_ * bin_cnt_y_);
  int x = lx_, y = ly_;
  int idxX = 0, idxY = 0;
  for (auto &bin : binStor_) {
    int sizeX = (x + bin_size_x_ > ux_) ? ux_ - x : bin_size_x_;
    int sizeY = (y + bin_size_y_ > uy_) ? uy_ - y : bin_size_y_;
    bin = Bin(idxX, idxY, x, y, x + sizeX, y + sizeY, target_density_);

    // move x, y coordinates.
    x += bin_size_x_;
    idxX += 1;
    if (x >= ux_) {
      y += bin_size_y_;
      x = lx_;
      idxY++;
      idxX = 0;
    }
    bins_.push_back(&bin);
  }


  // only initialized once
  updateBinsNonPlaceArea();
}
void Chip::NesterovPlacer::updateBinsNonPlaceArea() {
  for (auto &bin : bins_) {
    bin->setNonPlaceArea(0);
    bin->setNonPlaceAreaUnscaled(0);
  }

  for (auto &inst : nonPlaceInsts_) {
    std::pair<int, int> pairX = getMinMaxIdxX(inst);
    std::pair<int, int> pairY = getMinMaxIdxY(inst);
    for (int i = pairX.first; i < pairX.second; i++) {
      for (int j = pairY.first; j < pairY.second; j++) {
        Bin *bin = bins_[j * bin_cnt_x_ + i];

        // Note that nonPlaceArea should have scale-down with
        // target density.
        // See MS-replace paper
        //
        bin->addNonPlaceArea(getOverlapArea(bin, inst, db_database_->getChip()->getBlock()->getDbUnitsPerMicron())
                                 * bin->targetDensity());
        bin->addNonPlaceAreaUnscaled(getOverlapAreaUnscaled(bin, inst) * bin->targetDensity());
      }
    }
  }
}
float Chip::NesterovPlacer::getDensityCoordiLayoutInsideX(Instance *instance, float cx) {
  float adjVal = cx;  // adjusted value
  // will change base on each assigned binGrids.
  if (cx - instance->getDensityDeltaX() / 2 < die_pointer_->getLowerLeftX()) {
    adjVal = die_pointer_->getLowerLeftX() + instance->getDensityDeltaX() / 2;
  }
  if (cx + instance->getDensityDeltaX() / 2 > die_pointer_->getUpperRightX()) {
    adjVal = die_pointer_->getUpperRightX() - instance->getDensityDeltaX() / 2;
  }
  return adjVal;
}
float Chip::NesterovPlacer::getDensityCoordiLayoutInsideY(Instance *instance, float cy) {
  float adjVal = cy;  // adjusted value
  // will change base on each assigned binGrids.
  if (cy - instance->getDensityDeltaY() / 2 < die_pointer_->getLowerLeftY()) {
    adjVal = die_pointer_->getLowerLeftY() + instance->getDensityDeltaY() / 2;
  }
  if (cy + instance->getDensityDeltaY() / 2 > die_pointer_->getUpperRightY()) {
    adjVal = die_pointer_->getUpperRightY() - instance->getDensityDeltaY() / 2;
  }
  return adjVal;
}
void Chip::NesterovPlacer::updateGCellDensityCenterLocation(const vector<pair<float, float>> &coordinates) {
  for (int idx = 0; idx < coordinates.size(); ++idx) {
    pair<float, float> coordinate = coordinates.at(idx);
    Instance *instance = instance_pointers_.at(idx);
    instance->setDensityCenterLocation(coordinate.first, coordinate.second);
  }
  updateBinsCellDensityArea(instance_pointers_);
}
std::pair<int, int> Chip::NesterovPlacer::getMinMaxIdxX(Instance *inst) const {
  int lowerIdx = (inst->ly() - die_pointer_->getLowerLeftY()) / bin_size_y_;
  int upperIdx = (fastModulo((inst->uy() - lx()), bin_size_y_) == 0)
                 ? (inst->uy() - ly()) / bin_size_y_
                 : (inst->uy() - ly()) / bin_size_y_ + 1;
  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, bin_cnt_y_));
}
std::pair<int, int> Chip::NesterovPlacer::getMinMaxIdxY(Instance *inst) const {
  int lowerIdx = (inst->ly() - ly()) / bin_size_y_;
  int upperIdx = (fastModulo((inst->uy() - ly()), bin_size_y_) == 0)
                 ? (inst->uy() - ly()) / bin_size_y_
                 : (inst->uy() - ly()) / bin_size_y_ + 1;

  return std::make_pair(std::max(lowerIdx, 0), std::min(upperIdx, bin_cnt_y_));
}
int Chip::NesterovPlacer::fastModulo(const int input, const int ceil) { return input >= ceil ? input % ceil : input; }
int64_t Chip::NesterovPlacer::getOverlapArea(const Chip::NesterovPlacer::Bin *bin, Instance *inst, int dbu_per_micron) {
  int rectLx = max(bin->lx(), inst->lx()), rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), rectUy = min(bin->uy(), inst->uy());

  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  }

  if (inst->isMacro()) {
    float meanX = (inst->cx() - inst->lx()) / (float) dbu_per_micron;
    float meanY = (inst->cy() - inst->ly()) / (float) dbu_per_micron;

    // For the bivariate normal distribution, we are using
    // the shifted means of X and Y.
    // Sigma is used as the mean/4 for both dimensions
    const biNormalParameters i
        = {meanX,
           meanY,
           meanX / 4,
           meanY / 4,
           (rectLx - inst->lx()) / (float) dbu_per_micron,
           (rectLy - inst->ly()) / (float) dbu_per_micron,
           (rectUx - inst->lx()) / (float) dbu_per_micron,
           (rectUy - inst->ly()) / (float) dbu_per_micron};

    const float original = static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
    const float scaled = calculateBiVariateNormalCDF(i)
        * static_cast<float>(inst->ux() - inst->lx())
        * static_cast<float>(inst->uy() - inst->ly());

    // For heavily dense regions towards the center of the macro,
    // we are using an upper limit of 1.15*(overlap) between the macro
    // and the bin.
    if (scaled >= original) {
      return min<float>(scaled, original * 1.15);
    }
      // If the scaled value is smaller than the actual overlap
      // then use the original overlap value instead.
      // This is implemented to prevent cells from being placed
      // at the outer sides of the macro.
    else {
      return original;
    }
  } else {
    return static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
  }
}
float Chip::NesterovPlacer::calculateBiVariateNormalCDF(Chip::NesterovPlacer::biNormalParameters i) {
  const float x1 = (i.meanX - i.lx) / (sqrt(2) * i.sigmaX);
  const float x2 = (i.meanX - i.ux) / (sqrt(2) * i.sigmaX);

  const float y1 = (i.meanY - i.ly) / (sqrt(2) * i.sigmaY);
  const float y2 = (i.meanY - i.uy) / (sqrt(2) * i.sigmaY);

  return 0.25
      * (erf(x1) * erf(y1) + erf(x2) * erf(y2) - erf(x1) * erf(y2)
          - erf(x2) * erf(y1));
}
int64_t Chip::NesterovPlacer::getOverlapAreaUnscaled(const Chip::NesterovPlacer::Bin *bin, Instance *inst) {
  int rectLx = max(bin->lx(), inst->lx()), rectLy = max(bin->ly(), inst->ly()),
      rectUx = min(bin->ux(), inst->ux()), rectUy = min(bin->uy(), inst->uy());

  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  } else {
    return static_cast<int64_t>(rectUx - rectLx)
        * static_cast<int64_t>(rectUy - rectLy);
  }
}
void Chip::NesterovPlacer::updateDensitySize() {
  for (Instance *instance : instance_pointers_) {
    float scale_x = 0, scale_y = 0;
    float density_width = 0, density_height = 0;
    if (instance->dx() < REPLACE_SQRT2 * bin_size_x_) {
      scale_x = static_cast<float>(instance->dx()) / static_cast<float>(REPLACE_SQRT2 * bin_size_x_);
      density_width = REPLACE_SQRT2 * static_cast<float>(bin_size_x_);
    } else {
      scale_x = 1.0;
      density_width = instance->dx();
    }

    if (instance->dy() < REPLACE_SQRT2 * bin_size_y_) {
      scale_y = static_cast<float>(instance->dy())
          / static_cast<float>(REPLACE_SQRT2 * bin_size_y_);
      density_height = REPLACE_SQRT2 * static_cast<float>(bin_size_y_);
    } else {
      scale_y = 1.0;
      density_height = instance->dy();
    }

    instance->setDensitySize(density_width, density_height);
    instance->setDensityScale(scale_x * scale_y);
  }
}
void Chip::NesterovPlacer::updateDensityForceBin() {
  // copy density to utilize FFT
  for (auto &bin : bins_) {
    fft_->updateDensity(bin->x(), bin->y(), bin->density());
  }

  // do FFT
  fft_->doFFT();

  // update electroPhi and electroForce
  // update sum_phi_ for nesterov loop
  sum_phi_ = 0;
  for (auto &bin : bins_) {
    auto eForcePair = fft_->getElectroForce(bin->x(), bin->y());
    bin->setElectroForce(eForcePair.first, eForcePair.second);

    float electroPhi = fft_->getElectroPhi(bin->x(), bin->y());
    bin->setElectroPhi(electroPhi);

    sum_phi_ += electroPhi * static_cast<float>(bin->nonPlaceArea() + bin->instPlacedArea() + bin->fillerArea());
  }

}
void Chip::NesterovPlacer::updateWireLengthForceWA(double wlCoeffX, double wlCoeffY) {
  // clear all WA variables.
  for (Net *gNet : net_pointers_) {
    gNet->clearWaVars();
  }
  for (auto &gPin : pin_pointers_) {
    gPin->clearWaVars();
  }

  for (Net *&gNet : net_pointers_) {
    gNet->updateBox(die_pointer_->getDieId());
    vector<Pin *> pin_set;

    if (!gNet->isIntersected())
      pin_set = gNet->getConnectedPins();
    else {
      for (Pin *pin : gNet->getConnectedPins()) {
        if (pin->isInstancePin()) {
          if (pin->getInstance()->getDieId() == die_pointer_->getDieId())
            pin_set.push_back(pin);
        } else if (pin->isBlockPin()) {
          // TODO: more accurate method is needed.
          pin_set.push_back(pin);
        }
      }
    }

    for (Pin *pin : pin_set) {
      // The WA terms are shift invariant:
      //
      //   Sum(x_i * exp(x_i))    Sum(x_i * exp(x_i - C))
      //   -----------------    = -----------------
      //   Sum(exp(x_i))          Sum(exp(x_i - C))
      //
      // So we shift to keep the exponential from overflowing
      double expMinX = static_cast<double>(gNet->lx() - pin->cx()) * wlCoeffX;
      double expMaxX = static_cast<double>(pin->cx() - gNet->ux()) * wlCoeffX;
      double expMinY = static_cast<double>(gNet->ly() - pin->cy()) * wlCoeffY;
      double expMaxY = static_cast<double>(pin->cy() - gNet->uy()) * wlCoeffY;

      // min x
      if (expMinX > min_wire_length_force_bar_) {
        if (pin->isInstancePin())
          pin->setMinExpSumX(2 * fastExp(expMinX));
        else
          pin->setMinExpSumX(fastExp(expMinX));
        gNet->addWaExpMinSumX(pin->minExpSumX());
        gNet->addWaXExpMinSumX(pin->cx() * pin->minExpSumX());
      }

      // max x
      if (expMaxX > min_wire_length_force_bar_) {
        if (pin->isInstancePin())
          pin->setMaxExpSumX(2 * fastExp(expMaxX));
        else
          pin->setMaxExpSumX(fastExp(expMaxX));
        gNet->addWaExpMaxSumX(pin->maxExpSumX());
        gNet->addWaXExpMaxSumX(pin->cx() * pin->maxExpSumX());
      }

      // min y
      if (expMinY > min_wire_length_force_bar_) {
        if (pin->isInstancePin())
          pin->setMinExpSumY(2 * fastExp(expMinY));
        else
          pin->setMinExpSumY(fastExp(expMinY));
        gNet->addWaExpMinSumY(pin->minExpSumY());
        gNet->addWaYExpMinSumY(pin->cy() * pin->minExpSumY());
      }

      // max y
      if (expMaxY > min_wire_length_force_bar_) {
        if (pin->isInstancePin())
          pin->setMaxExpSumY(2 * fastExp(expMaxY));
        else
          pin->setMaxExpSumY(fastExp(expMaxY));
        gNet->addWaExpMaxSumY(pin->maxExpSumY());
        gNet->addWaYExpMaxSumY(pin->cy() * pin->maxExpSumY());
      }
    }
  }

}
float Chip::NesterovPlacer::nesterovInstsArea() const {
  return std_instances_area_ + static_cast<int64_t>(round(macro_instances_area_ * target_density_));
}
void Chip::NesterovPlacer::updateWireLengthCoef(float overflow) {
  if (overflow > 1.0) {
    wire_length_coefficient_x_ = wire_length_coefficient_y_ = 0.1;
  } else if (overflow < 0.1) {
    wire_length_coefficient_x_ = wire_length_coefficient_y_ = 10.0;
  } else {
    wire_length_coefficient_x_ = wire_length_coefficient_y_
        = 1.0 / pow(10.0, (overflow - 0.1) * 20 / 9.0 - 1.0);
  }
  wire_length_coefficient_x_ *= base_wire_length_coefficient_;
  wire_length_coefficient_y_ *= base_wire_length_coefficient_;
}
float Chip::NesterovPlacer::getPhiCoef(float scaledDiffHpwl) {
  float retCoef = (scaledDiffHpwl < 0)
                  ? max_phi_coef_
                  : max_phi_coef_
                      * pow(max_phi_coef_, scaledDiffHpwl * -1.0);
  retCoef = std::max(minPhiCoef, retCoef);
  return retCoef;
}
int64_t Chip::NesterovPlacer::getHpwl() {
  int64_t hpwl = 0;
  for (auto &gNet : net_pointers_) {
    gNet->updateBox(this->die_pointer_->getDieId());
    hpwl += gNet->hpwl();
  }
  return hpwl;
}
void Chip::NesterovPlacer::updateNextIter() {
  // swap vector pointers
  std::swap(prev_slp_coordinates_, cur_slp_coordinates_);
  std::swap(prev_slp_wire_length_grads_, cur_slp_wire_length_grads_);
  std::swap(prev_slp_density_grads_, cur_slp_density_grads_);
  std::swap(prev_slp_sum_grads_, cur_slp_sum_grads_);

  // Prevent locked instances from moving
  const auto &gCells = instance_pointers_;
  for (size_t k = 0; k < gCells.size(); ++k) {
    if (gCells[k]->isInstance() && gCells[k]->isLocked()) {
      next_slp_coordinates_[k] = cur_slp_coordinates_[k];
      next_slp_wire_length_grads_[k] = cur_slp_wire_length_grads_[k];
      next_slp_density_grads_[k] = cur_slp_density_grads_[k];
      next_slp_sum_grads_[k] = cur_slp_sum_grads_[k];

      next_coordinates_[k] = cur_coordinates_[k];
    }
  }

  std::swap(cur_slp_coordinates_, next_slp_coordinates_);
  std::swap(cur_slp_wire_length_grads_, next_slp_wire_length_grads_);
  std::swap(cur_slp_density_grads_, next_slp_density_grads_);
  std::swap(cur_slp_sum_grads_, next_slp_sum_grads_);

  std::swap(cur_coordinates_, next_coordinates_);

  sum_overflow_ = static_cast<float>(overflow_area_) / static_cast<float>(nesterovInstsArea());

  sum_overflow_unscaled_ = static_cast<float>(overflow_area_unscaled_) / static_cast<float>(nesterovInstsArea());

  updateWireLengthCoef(sum_overflow_);
  int64_t hpwl = getHpwl();

  float phiCoef = getPhiCoef(static_cast<float>(hpwl - prev_hpwl_) / referenceHpwl);

  prev_hpwl_ = hpwl;
  density_penalty_ *= phiCoef;

/*
      // for routability densityPenalty recovery
      if (rb_->numCall() == 0) {
        density_penalty_storage_.push_back(density_penalty_);
      }
*/
}
pair<float, float> Chip::NesterovPlacer::getWireLengthPreconditioner(Instance *instance) {
  // original function: getWireLengthPreconditioner
  int binding_nums = 0;
  for (Pin *pin : instance->getPins()) {
    if (pin->getNet())
      binding_nums += 1;
  }
  return pair<float, float>{static_cast<float>(binding_nums), static_cast<float>(binding_nums)};
}
pair<float, float> Chip::NesterovPlacer::getDensityPreconditioner(Instance *gCell) {
  float areaVal
      = static_cast<float>(gCell->dx()) * static_cast<float>(gCell->dy());

  return pair<float, float>{areaVal, areaVal};

}
std::pair<int, int> Chip::NesterovPlacer::getDensityMinMaxIdxX(Instance *gcell) {
  int lowerIdx = (gcell->dLx() - lx()) / bin_size_x_;
  int upperIdx = (fastModulo((gcell->dUx() - lx()), bin_size_x_) == 0)
                 ? (gcell->dUx() - lx()) / bin_size_x_
                 : (gcell->dUx() - lx()) / bin_size_x_ + 1;

  upperIdx = std::min(upperIdx, bin_cnt_x_);
  return std::make_pair(lowerIdx, upperIdx);
}
std::pair<int, int> Chip::NesterovPlacer::getDensityMinMaxIdxY(Instance *gcell) {
  int lowerIdx = (gcell->dLy() - ly()) / bin_size_y_;
  int upperIdx = (fastModulo((gcell->dUy() - ly()), bin_size_y_) == 0)
                 ? (gcell->dUy() - ly()) / bin_size_y_
                 : (gcell->dUy() - ly()) / bin_size_y_ + 1;

  upperIdx = std::min(upperIdx, bin_cnt_y_);
  return std::make_pair(lowerIdx, upperIdx);
}
float Chip::NesterovPlacer::getOverlapDensityArea(Chip::NesterovPlacer::Bin *bin, Instance *cell) {
  int rectLx = max(bin->lx(), static_cast<int>(cell->dLx()));
  int rectLy = max(bin->ly(), static_cast<int>(cell->dLy()));
  int rectUx = min(bin->ux(), static_cast<int>(cell->dUx()));
  int rectUy = min(bin->uy(), static_cast<int>(cell->dUy()));
  if (rectLx >= rectUx || rectLy >= rectUy) {
    return 0;
  } else {
    return static_cast<float>(rectUx - rectLx)
        * static_cast<float>(rectUy - rectLy);
  }
}
pair<float, float> Chip::NesterovPlacer::getDensityGradient(Instance *gCell) {
  std::pair<int, int> pairX = getDensityMinMaxIdxX(gCell);
  std::pair<int, int> pairY = getDensityMinMaxIdxY(gCell);

  pair<float, float> electroForce;

  for (int i = pairX.first; i < pairX.second; i++) {
    for (int j = pairY.first; j < pairY.second; j++) {
      Bin *bin = bins_.at(j * bin_cnt_x_ + i);
      float overlapArea
          = getOverlapDensityArea(bin, gCell) * gCell->densityScale();

      electroForce.first += overlapArea * bin->electroForceX();
      electroForce.second += overlapArea * bin->electroForceY();
    }
  }
  return electroForce;

}
void Chip::NesterovPlacer::updateGradients(vector<pair<float, float>> &sumGrads,
                                           vector<pair<float, float>> &wireLengthGrads,
                                           vector<pair<float, float>> &densityGrads) {
  wire_length_grad_sum_ = 0;
  density_grad_sum_ = 0;

  float gradSum = 0;

  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *cell = instance_pointers_.at(i);
    wireLengthGrads[i] = getWireLengthGradientWA(cell, wire_length_coefficient_x_, wire_length_coefficient_y_);
    densityGrads[i] = getDensityGradient(cell);

    // Different compiler has different results on the following formula.
    // e.g. wire_length_grad_sum_ += fabs(~~.x) + fabs(~~.y);
    //
    // To prevent instability problem,
    // I partitioned the fabs(~~.x) + fabs(~~.y) as two terms.
    //
    wire_length_grad_sum_ += fabs(wireLengthGrads[i].first);
    wire_length_grad_sum_ += fabs(wireLengthGrads[i].second);

    density_grad_sum_ += fabs(densityGrads[i].first);
    density_grad_sum_ += fabs(densityGrads[i].second);

    sumGrads[i].first = wireLengthGrads[i].first + density_penalty_ * densityGrads[i].first;
    sumGrads[i].second = wireLengthGrads[i].second + density_penalty_ * densityGrads[i].second;

    pair<float, float> wire_length_preconditioner = getWireLengthPreconditioner(cell);
    pair<float, float> density_preconditioner = getDensityPreconditioner(cell);

    pair<float, float> sum_preconditioner(
        wire_length_preconditioner.first + density_penalty_ * density_preconditioner.first,
        wire_length_preconditioner.second + density_penalty_ * density_preconditioner.second);

    if (sum_preconditioner.first <= minPreconditioner) {
      sum_preconditioner.first = minPreconditioner;
    }

    if (sum_preconditioner.second <= minPreconditioner) {
      sum_preconditioner.second = minPreconditioner;
    }

    sumGrads[i].first /= sum_preconditioner.first;
    sumGrads[i].second /= sum_preconditioner.second;

    gradSum += fabs(sumGrads[i].first) + fabs(sumGrads[i].second);
  }

  // sometimes wirelength gradient is zero when design is too small
  if (wire_length_grad_sum_ == 0
      && recursion_cnt_wl_coef_ < maxRecursionInitSLPCoef) {
    wire_length_coefficient_x_ *= 0.5;
    wire_length_coefficient_y_ *= 0.5;
    base_wire_length_coefficient_ *= 0.5;

    // update WL forces
    updateWireLengthForceWA(wire_length_coefficient_x_, wire_length_coefficient_y_);

    // recursive call again with smaller wirelength coef
    recursion_cnt_wl_coef_++;
  }

  // divergence detection on
  // Wirelength / density gradient calculation
  if (isnan(wire_length_grad_sum_) || isinf(wire_length_grad_sum_)
      || isnan(density_grad_sum_) || isinf(density_grad_sum_)) {
    is_diverged_ = true;
    diverge_msg_ = "RePlAce diverged at wire/density gradient Sum.";
    diverge_code_ = 306;
  }
}
float Chip::NesterovPlacer::getStepLength() {
  float coordiDistance = getDistance(prev_slp_coordinates_, cur_slp_coordinates_);
  float gradDistance = getDistance(prev_slp_sum_grads_, cur_slp_sum_grads_);
  return coordiDistance / gradDistance;
}
float Chip::NesterovPlacer::getDistance(vector<pair<float, float>> a, vector<pair<float, float>> b) {
  float sumDistance = 0.0f;
  for (size_t i = 0; i < a.size(); i++) {
    sumDistance += (a[i].first - b[i].first) * (a[i].first - b[i].first);
    sumDistance += (a[i].second - b[i].second) * (a[i].second - b[i].second);
  }

  return sqrt(sumDistance / (2.0 * a.size()));
}
pair<float, float> Chip::NesterovPlacer::getWireLengthGradientWA(Instance *gCell,
                                                                 float wlCoeffX,
                                                                 float wlCoeffY) const {
  pair<float, float> gradientPair;

  for (auto &gPin : gCell->getPins()) {
    if (gPin->getNet() == nullptr)
      // pass the floating pins
      continue;

    auto tmpPair = getWireLengthGradientPinWA(gPin, wlCoeffX, wlCoeffY);

    // apply timing/custom net weight
    tmpPair.first *= gPin->getNet()->totalWeight();
    tmpPair.second *= gPin->getNet()->totalWeight();

    gradientPair.first += tmpPair.first;
    gradientPair.second += tmpPair.second;
  }

  // return sum
  assert(isnan(gradientPair.first) == false);
  assert(isnan(gradientPair.second) == false);
  return gradientPair;
}
pair<float, float> Chip::NesterovPlacer::getWireLengthGradientPinWA(Pin *gPin, float wlCoeffX, float wlCoeffY) {
  double gradientMinX = 0, gradientMinY = 0;
  double gradientMaxX = 0, gradientMaxY = 0;

  // min x
  if (gPin->hasMinExpSumX()) {
    // from Net.
    double waExpMinSumX = gPin->getNet()->waExpMinSumX();
    double waXExpMinSumX = gPin->getNet()->waXExpMinSumX();

    gradientMinX
        = static_cast<double>(waExpMinSumX * (gPin->minExpSumX() * (1.0 - wlCoeffX * static_cast<double>(gPin->cx())))
        + wlCoeffX * gPin->minExpSumX() * waXExpMinSumX) / (waExpMinSumX * waExpMinSumX);
  }

  // max x
  if (gPin->hasMaxExpSumX()) {
    double waExpMaxSumX = gPin->getNet()->waExpMaxSumX();
    double waXExpMaxSumX = gPin->getNet()->waXExpMaxSumX();

    gradientMaxX
        = static_cast<double>(waExpMaxSumX * (gPin->maxExpSumX() * (1.0 + wlCoeffX * static_cast<double>(gPin->cx())))
        - wlCoeffX * gPin->maxExpSumX() * waXExpMaxSumX) / (waExpMaxSumX * waExpMaxSumX);
  }

  // min y
  if (gPin->hasMinExpSumY()) {
    double waExpMinSumY = gPin->getNet()->waExpMinSumY();
    double waYExpMinSumY = gPin->getNet()->waYExpMinSumY();

    gradientMinY
        = static_cast<double>(waExpMinSumY * (gPin->minExpSumY() * (1.0 - wlCoeffY * static_cast<double>(gPin->cy())))
        + wlCoeffY * gPin->minExpSumY() * waYExpMinSumY) / (waExpMinSumY * waExpMinSumY);
  }

  // max y
  if (gPin->hasMaxExpSumY()) {
    double waExpMaxSumY = gPin->getNet()->waExpMaxSumY();
    double waYExpMaxSumY = gPin->getNet()->waYExpMaxSumY();

    gradientMaxY
        = static_cast<double>(waExpMaxSumY * (gPin->maxExpSumY() * (1.0 + wlCoeffY * static_cast<double>(gPin->cy())))
        - wlCoeffY * gPin->maxExpSumY() * waYExpMaxSumY) / (waExpMaxSumY * waExpMaxSumY);
  }

  assert(!isnan(gradientMaxX));
  assert(!isnan(gradientMaxY));
  assert(!isnan(gradientMinX));
  assert(!isnan(gradientMinY));
  return pair<float, float>{gradientMinX - gradientMaxX, gradientMinY - gradientMaxY};
}
void Chip::NesterovPlacer::initSLPStepsVars() {
  const int instance_num = instance_pointers_.size();
  cur_slp_coordinates_.resize(instance_num);
  cur_slp_wire_length_grads_.resize(instance_num);
  cur_slp_density_grads_.resize(instance_num);
  cur_slp_sum_grads_.resize(instance_num);
  next_slp_coordinates_.resize(instance_num);
  next_slp_wire_length_grads_.resize(instance_num);
  next_slp_density_grads_.resize(instance_num);
  next_slp_sum_grads_.resize(instance_num);
  prev_slp_coordinates_.resize(instance_num);
  prev_slp_wire_length_grads_.resize(instance_num);
  prev_slp_density_grads_.resize(instance_num);
  prev_slp_sum_grads_.resize(instance_num);
  cur_coordinates_.resize(instance_num);
  next_coordinates_.resize(instance_num);
  init_coordinates_.resize(instance_num);
}
void Chip::NesterovPlacer::updateBinsCellDensityArea(vector<Instance *> cells) {
  // clear the Bin-area info
  for (auto &bin : bins_) {
    bin->setInstPlacedArea(0);
    bin->setInstPlacedAreaUnscaled(0);
    bin->setFillerArea(0);
  }

  for (auto &cell : cells) {
    std::pair<int, int> pairX = getDensityMinMaxIdxX(cell);
    std::pair<int, int> pairY = getDensityMinMaxIdxY(cell);

    // The following function is critical runtime hotspot
    // for global placer.
    //
    if (cell->isInstance()) {
      // macro should have
      // scale-down with target-density
      if (cell->isMacroInstance()) {
        for (int i = pairX.first; i < pairX.second; i++) {
          for (int j = pairY.first; j < pairY.second; j++) {
            Bin *bin = bins_[j * bin_cnt_x_ + i];

            const float scaledAvea = getOverlapDensityArea(bin, cell)
                * cell->densityScale()
                * bin->targetDensity();
            bin->addInstPlacedArea(scaledAvea);
            bin->addInstPlacedAreaUnscaled(scaledAvea);
          }
        }
      }
        // normal cells
      else if (cell->isStdInstance()) {
        for (int i = pairX.first; i < pairX.second; i++) {
          for (int j = pairY.first; j < pairY.second; j++) {
            Bin *bin = bins_[j * bin_cnt_x_ + i];
            const float scaledArea
                = getOverlapDensityArea(bin, cell) * cell->densityScale();
            bin->addInstPlacedArea(scaledArea);
            bin->addInstPlacedAreaUnscaled(scaledArea);
          }
        }
      }
    } else if (cell->isFiller()) {
      for (int i = pairX.first; i < pairX.second; i++) {
        for (int j = pairY.first; j < pairY.second; j++) {
          Bin *bin = bins_[j * bin_cnt_x_ + i];
          bin->addFillerArea(getOverlapDensityArea(bin, cell)
                                 * cell->densityScale());
        }
      }
    }
  }

  overflow_area_ = 0;
  overflow_area_unscaled_ = 0;
  // update density and overflowArea
  // for nesterov use and FFT library
  for (auto &bin : bins_) {
    int64_t binArea = bin->binArea();
    const float scaledBinArea = static_cast<float>(binArea * bin->targetDensity());
    bin->setDensity((static_cast<float>(bin->instPlacedArea())
        + static_cast<float>(bin->fillerArea()) + static_cast<float>(bin->nonPlaceArea())) / scaledBinArea);

    overflow_area_ += std::max(0.0f,
                               static_cast<float>(bin->instPlacedArea()) + static_cast<float>(bin->nonPlaceArea())
                                   - scaledBinArea);

    overflow_area_unscaled_ += std::max(
        0.0f,
        static_cast<float>(bin->instPlacedAreaUnscaled())
            + static_cast<float>(bin->nonPlaceAreaUnscaled()) - scaledBinArea);
  }
}
void Chip::NesterovPlacer::setDensityValuesAsDefault() {
  for (Instance *instance : instance_pointers_) {
    instance->setDensityValueAsDefault();
  }
}
void Chip::NesterovPlacer::updateDensityCoordiLayoutInside(Instance *gCell) {
  float targetLx = gCell->dLx();
  float targetLy = gCell->dLy();

  if (targetLx < lx()) {
    targetLx = lx();
  }

  if (targetLy < ly()) {
    targetLy = ly();
  }

  if (targetLx + gCell->getDensityDeltaX() > ux()) {
    targetLx = ux() - gCell->getDensityDeltaX();
  }

  if (targetLy + gCell->getDensityDeltaY() > uy()) {
    targetLy = uy() - gCell->getDensityDeltaY();
  }
  gCell->setDensityLocation(targetLx, targetLy);

}
void Chip::NesterovPlacer::updateInitialPrevSLPCoordi() {
  for (size_t i = 0; i < instance_pointers_.size(); i++) {
    Instance *curGCell = instance_pointers_.at(i);

    float prevCoordiX
        = cur_slp_coordinates_[i].first
            - initialPrevCoordiUpdateCoef * cur_slp_sum_grads_[i].first;

    float prevCoordiY
        = cur_slp_coordinates_[i].second
            - initialPrevCoordiUpdateCoef * cur_slp_sum_grads_[i].second;

    pair<float, float> newCoordi(getDensityCoordiLayoutInsideX(curGCell, prevCoordiX),
                                 getDensityCoordiLayoutInsideY(curGCell, prevCoordiY));

    prev_slp_coordinates_[i] = newCoordi;
  }

}
void Chip::NesterovPlacer::updateDB() {
  for (Instance *instance : instance_pointers_) {
    instance->setCoordinate(instance->getDensityCenterX(), instance->getDensityCenterY());
  }
}
void Chip::NesterovPlacer::setDebugMode(bool debug_mode) {
  debug_mode_ = debug_mode;
}
void Chip::NesterovPlacer::drawDie(const string &filename) {
  assert(debug_mode_);
  int fixed_die_height = 500;
  int scale_factor = die_pointer_->getHeight() / fixed_die_height;

  uint die_w = die_pointer_->getWidth() / scale_factor;
  uint die_h = die_pointer_->getHeight() / scale_factor;
  Drawer drawer(die_w, die_h);

  // Draw cells and fillers
  for (int i = 0; i < instance_pointers_.size(); ++i) {
    Instance *instance = instance_pointers_.at(i);
    pair<float, float> coordinate = cur_coordinates_.at(i);
    int instance_width = instance->dx();
    int instance_height = instance->dy();

    int ll_x = static_cast<int>(coordinate.first / scale_factor);
    int ll_y = static_cast<int>(coordinate.second / scale_factor);
    int ur_x = ll_x + static_cast<int>(instance_width / scale_factor);
    int ur_y = ll_y + static_cast<int>(instance_height / scale_factor);
    if (!instance->isFiller())
      drawer.drawCell(ll_x, ll_y, ur_x, ur_y);
    else
      drawer.drawFiller(ll_x, ll_y, ur_x, ur_y);
  }
  drawer.saveImg(filename);
}
double fastExp(float a) {
  a = 1.0f + a / 1024.0f;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  a *= a;
  return a;
}
Chip::NesterovPlacer::Drawer::Drawer(uint width, uint height) {
  width_ = width;
  height_ = height;
  image_ = new Image(width_, height_, 1, 3, 255);
}
Chip::NesterovPlacer::Drawer::~Drawer() {
  delete image_;
}
void Chip::NesterovPlacer::Drawer::setCellColor(const unsigned char *cell_color) {
  cell_color_ = cell_color;
}
void Chip::NesterovPlacer::Drawer::setFillerColor(const unsigned char *filler_color) {
  filler_color_ = filler_color;
}
void Chip::NesterovPlacer::Drawer::drawCell(int ll_x, int ll_y, int ur_x, int ur_y) {
  if (ll_x == ur_x)
    ur_x += 1;
  if (ll_y == ur_y)
    ur_y += 1;
  image_->draw_rectangle(ll_x, ll_y, ur_x, ur_y, Color::DIM_GRAY);
  image_->draw_rectangle(ll_x + 1, ll_y + 1, ur_x - 1, ur_y - 1, cell_color_);
}
void Chip::NesterovPlacer::Drawer::drawFiller(int ll_x, int ll_y, int ur_x, int ur_y) {
  if (ll_x == ur_x)
    ur_x += 1;
  if (ll_y == ur_y)
    ur_y += 1;
  image_->draw_rectangle(ll_x, ll_y, ur_x, ur_y, Color::DIM_GRAY);
  image_->draw_rectangle(ll_x + 1, ll_y + 1, ur_x - 1, ur_y - 1, filler_color_);
}
void Chip::NesterovPlacer::Drawer::saveImg(const string &file_name) {
  string save_file_name = file_path_ + file_name + ".bmp";
  image_->save(save_file_name.c_str());
}

// etc //
__attribute__((unused)) void Chip::dbTutorial() const {
  cout << this->db_database_->getChip()->getBlock()->getBBox()->getDX() << endl;

  dbBlock *block = db_database_->getChip()->getBlock();
  for (int i = 0; i < 4; ++i) {
    cout << endl;
  }
  // instance cloning into circuit variables
  dbSet<dbInst> db_instances = block->getInsts();
  int cellIdx = 0;

  for (dbInst *db_instance : db_instances) {
    if (cellIdx % 10 == 0) {
      cout << "db_instance->getMaster()->getLib()->getDbUnitsPerMicron(): "
           << db_instance->getMaster()->getLib()->getDbUnitsPerMicron() << endl;
      int x, y;
      auto type = db_instance->getMaster()->getType();
      if (!type.isCore() && !type.isBlock()) continue;

      cout << "db_instance->getName(): " << db_instance->getName() << endl;
      db_instance->getOrigin(x, y);
      cout << "getOrigin: " << x << " " << y << endl;
      cout << "db_instance->getMaster()->getName(): " << db_instance->getMaster()->getName() << endl;
      cout << "db_instance->getMaster()->getHeight(): " << db_instance->getMaster()->getHeight() << endl;
      cout << "db_instance->getMaster()->getWidth(): " << db_instance->getMaster()->getWidth() << endl;
      cout << "db_instance->getMaster()->getLib()->getName(): " << db_instance->getMaster()->getLib()->getName()
           << endl << endl;

      for (auto masterTerminal : db_instance->getMaster()->getMTerms()) {
        cout << endl;
        cout << "terminal name: " << masterTerminal->getName() << endl;
        cout << "masterTerminal->getSigType().getString(): " << masterTerminal->getSigType().getString() << endl;
        cout << "masterTerminal->getIoType().getString(): " << masterTerminal->getIoType().getString() << endl;
        cout << "masterTerminal->getIndex(): " << masterTerminal->getIndex() << endl;
        cout << "for (auto pin: masterTerminal->getMPins())" << endl;
        for (auto pin : masterTerminal->getMPins()) {
          for (auto box : pin->getGeometry()) {
            cout << "\t" << box->xMin() << " " << box->xMax() << "  " << box->yMin() << " " << box->yMax() << endl;
          }
        }
      }
      cout << endl;

      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        if (instanceTerminal->getNet())
          cout << instanceTerminal->getNet()->getName() << endl;
        else
          cout << "None net." << endl;
        cout << " " << instanceTerminal->getMTerm()->getName()
             << " " << instanceTerminal->getIoType().getString() << " " << x << " " << y << endl;;
      }

      cout << endl << endl;
      cout << "get origin: " << x << " " << y << endl;
      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        cout << x << " " << y << endl;
      }
      db_instance->setOrigin(100, 200);
      db_instance->getOrigin(x, y);
      cout << "get origin: " << x << " " << y << endl;
      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        cout << x << " " << y << endl;
      }

      cout << "instance end." << endl << endl;
      cellIdx++;
    }
  }
/*

  cout << endl << endl << endl << endl;
  cout << "OpenDB Tutorial starts." << endl;

  dbBlock *block = parser_.db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();

  // How to access all instances in database, and print the instance name
  int numOfInstance = 0;
  for (dbInst *db_instance : db_instances) {
    if (!db_instance->getName().empty())
      numOfInstance++;
  }
  cout << "The number of Instances: " << numOfInstance << endl;

  // Let's access only the first one instance, and explore some other methods.
  dbInst *db_inst = *db_instances.begin();
  cout << "The name of the instance: " << db_inst->getName() << endl;
  // As the library of the instance, we can know the library name, and cell width and height.
  dbMaster *lib = db_inst->getMaster();
  cout << "The name of the library of the instance: " << lib->getName() << endl;
  cout << "The width of the instance: " << db_inst->getMaster()->getWidth() << endl;
  cout << "The height of the instance: " << db_inst->getMaster()->getHeight() << endl << endl;

  cout << "The pin libraries of the instance" << endl;
  for (auto library_terminals : db_inst->getMaster()->getMTerms()) {
    cout << " ";
    cout << "Name: " << library_terminals->getName() << "\t";
    cout << "Signal Type: " << library_terminals->getSigType().getString() << "\t";
    cout << "Terminal Index: " << library_terminals->getIndex() << "\t";
    cout << endl;
  }
  cout << "The number of the terminals of instance: " << db_inst->getMaster()->getMTermCount() << endl;
  cout << endl;

  cout << "The pin data" << endl;
  for (auto instance_terminal : db_inst->getITerms()) {
    cout << " ";
    cout << "Terminal Name: " << instance_terminal->getMTerm()->getName() << "\t";
    cout << "Signal Type: " << instance_terminal->getSigType().getString() << "\t";
    cout << "Terminal Index: " << instance_terminal->getMTerm()->getIndex() << endl;
    cout << "  Connected Net: ";
    if (instance_terminal->getNet() != nullptr) {
      cout << instance_terminal->getNet()->getName() << endl;
      if (instance_terminal->getNet()->getName() != "clk") {
        cout << "   connected Pins correspond to the net" << endl;
        for (auto terminals_in_net : instance_terminal->getNet()->getITerms()) {
          cout << "   "
               << terminals_in_net->getMTerm()->getName() << " in "
               << terminals_in_net->getInst()->getName() << "("
               << terminals_in_net->getInst()->getMaster()->getName() << ") instance" << endl;
        }
      }
    } else {
      cout << "None" << endl;
    }
  }

  // # 2. traverse nets
  for (int i = 0; i < 3; ++i) {
    cout << endl;
  }
  // access the nets in the circuit
  cout << "Signal Type\tIoType\tTerminal Name" << endl;
  for (auto net : block->getNets()) {
    for (auto terminal : net->getITerms()) {
    }
    for (auto terminal : net->getBTerms()) {
      cout << terminal->getSigType().getString() << "\t";
      cout << terminal->getIoType().getString() << "\t";
      cout << terminal->getName() << endl;
      for (auto bPin : terminal->getBPins()) {
        cout << " bPin ID: " << bPin->getId() << endl;
        for (auto box : bPin->getBoxes()) {
          cout << "  box xMin: " << box->xMin() << endl;
          cout << "  box xMax: " << box->xMax() << endl;
          cout << "  box DX: " << box->getDX() << endl;
          cout << "  box yMin: " << box->yMin() << endl;
          cout << "  box yMax: " << box->yMax() << endl;
          cout << "  box DY: " << box->getDY() << endl;
          cout << "  box width: " << box->getWidth() << endl;
          cout << "  box dir: " << box->getDir() << endl;
          if (box->getBlockVia())*/
/* meaningless in this case*//*

            cout << "  via name: " << box->getBlockVia()->getName() << endl;
        }

      }
    }
  }


  // # 3. Get Die information
  odb::Rect Die, Core;
  block->getDieArea(Die);
  block->getCoreArea(Core);
  cout << Die.dx() << " " << Die.dy() << endl;
  cout << Core.dx() << " " << Core.dy() << endl;
*/

}
} // VLSI_backend