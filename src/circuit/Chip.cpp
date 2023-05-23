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
  parser_.setLoggerPtr(&logger_);
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
ulong Chip::getHPWL() {
  ulong HPWL = 0;
  for (Net *net : net_pointers_) {
    net->updateBox();
    HPWL += net->getHPWL();
  }
  return HPWL;
}
int Chip::getUnitOfMicro() const {
  return db_database_->getTech()->getDbUnitsPerMicron();
}
void Chip::drawDies(const string &pseudo_die_name,
                    const string &top_die_name,
                    const string &bottom_die_name,
                    int scale_factor,
                    bool as_dot,
                    bool draw_same_canvas) {
  if (instance_pointers_.size() > 10e6) {
    scale_factor = 1000;
  } else if (instance_pointers_.size() > 10e3) {
    scale_factor = 100;
  } else {
    scale_factor = 10;
  }

  if (!data_storage_.hybrid_bonds.empty()) {
    uint top_die_w = die_pointers_.at(0)->getWidth() / scale_factor;
    uint top_die_h = die_pointers_.at(0)->getHeight() / scale_factor;
    uint bottom_die_w = die_pointers_.at(1)->getWidth() / scale_factor;
    uint bottom_die_h = die_pointers_.at(1)->getHeight() / scale_factor;

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
    top_die.setDieId(1);
    bottom_die.setDieId(2);

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
      top_die.saveImg(pseudo_die_name);
    else {
      top_die.saveImg(top_die_name);
      bottom_die.saveImg(bottom_die_name);
    }
  } else {
    uint pseudo_die_w = die_pointers_.at(0)->getWidth() / scale_factor;
    uint pseudo_die_h = die_pointers_.at(0)->getHeight() / scale_factor;

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

    pseudo_die.saveImg(pseudo_die_name);
  }

}
} // VLSI_backend