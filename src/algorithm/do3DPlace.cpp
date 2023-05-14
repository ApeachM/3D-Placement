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
#include "NesterovPlacer.h"
#include "InitialPlacer.h"
#include <random>
#include <algorithm>
namespace VLSI_backend {
void Chip::do3DPlace() {

  // 0. target density setting
  /* manually setting in code level */
  vector<double> densities;
  densities.push_back(2.0); // pseudo die util = top die util + bottom die util
  densities.push_back(1.0);
  densities.push_back(1.0);
  setTargetDensity(densities);

  // 1. do3DPlace the cells in the pseudo die
  this->normalPlacement();

  // 2. partition
  this->partition();  // temporary code

  // 3. hybrid bond generate and placement
  this->generateHybridBonds();

  // 4. placement synchronously
  this->placement2DieSynchronously();
}

void Chip::normalPlacement() {
  doInitialPlace();
  doNestrovPlace();
  this->drawDies();
}

void Chip::partition() {
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
  data_storage_.hybrid_bond_pins.reserve(net_pointers_.size());

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
      Instance hybrid_bond_object;
      Pin hybrid_bond_pin_object;
      hybrid_bond_object.setAsHybridBond();
      hybrid_bond_pin_object.setAsHybridBondPin();

      // store them in storage
      data_storage_.hybrid_bonds.push_back(hybrid_bond_object);
      data_storage_.hybrid_bond_pins.push_back(hybrid_bond_pin_object);

      Instance *hybrid_bond = &data_storage_.hybrid_bonds.at(data_storage_.hybrid_bonds.size() - 1);
      Pin *hybrid_bond_pin = &data_storage_.hybrid_bond_pins.at(data_storage_.hybrid_bond_pins.size() - 1);

      // set name
      hybrid_bond->setInstName("hybrid_bond_" + to_string(hybrid_num));

      // link them
      hybrid_bond->setHybridBondPin(hybrid_bond_pin); // pin and instance
      hybrid_bond_pin->setHybridBond(hybrid_bond);
      net->setHybridBondPin(hybrid_bond_pin); // pin and net
      hybrid_bond_pin->setIntersectedNet(net);

      // set coordinate of hybrid bond
      // p.s. net box would be updated in the first placement phase (Nestrov in virtual die)
      int center_x = static_cast<int>((net->ux() + net->lx()) / 2);
      int center_y = static_cast<int>((net->uy() + net->ly()) / 2);
      hybrid_bond->setCoordinate(center_x, center_y);


      // for iteration
      instance_pointers_.push_back(hybrid_bond);
      pin_pointers_.push_back(hybrid_bond_pin);
    }
  }

  assert(hybrid_num == data_storage_.hybrid_bonds.size());
  cout << "Hybrid bond #: " << data_storage_.hybrid_bonds.size() << endl;
  this->drawDies("after_hybrid_bond_generation");
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
    if (die_id == -1) {
      // if die_id == -1, it is a hybrid bond
      assert(instance->isHybridBond());
      dieVar1.instance_pointers.push_back(instance);
      dieVar2.instance_pointers.push_back(instance);
    } else if (die_id == 1) {
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
    if (die_id == -1) {
      // intersected
      dieVar1.net_pointers_.push_back(net);
      dieVar2.net_pointers_.push_back(net);
    } else if (die_id == 1)
      // the net is on die 1
      dieVar1.net_pointers_.push_back(net);
    else if (die_id == 2)
      // the net is on die 2
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

    if (die_id == -1) {
      // intersected
      dieVar1.pin_pointers_.push_back(pin);
      dieVar2.pin_pointers_.push_back(pin);
    } else if (die_id == 1) {
      dieVar1.pin_pointers_.push_back(pin);
    } else if (die_id == 2) {
      dieVar2.pin_pointers_.push_back(pin);
    } else
      assert(0);
  }

  // skip the pad pointers because the nestrov optimizer doesn't use them.


  /////////////////////////////////////////////////////////////////////////
  NesterovPlacer nestrov_placer1(
      this->db_database_,
      dieVar1.instance_pointers,
      dieVar1.net_pointers_,
      dieVar1.pin_pointers_,
      dieVar1.pad_pointers_,
      this->die_pointers_.at(1)
  );
  NesterovPlacer nestrov_placer2(
      this->db_database_,
      dieVar2.instance_pointers,
      dieVar2.net_pointers_,
      dieVar2.pin_pointers_,
      dieVar2.pad_pointers_,
      this->die_pointers_.at(2)
  );

  nestrov_placer1.initNestrovPlace();
  nestrov_placer2.initNestrovPlace();

  if (nestrov_placer1.getMaxNesterovIter() != nestrov_placer2.getMaxNesterovIter())
    assert(0);

  for (int i = 0; i < nestrov_placer1.getMaxNesterovIter(); ++i) {
    int nestrov_iter1, nestrov_iter2;
    nestrov_iter1 = nestrov_placer1.doNestrovPlace(i, true);
    nestrov_iter2 = nestrov_placer2.doNestrovPlace(i, true);
    if (nestrov_iter1 >= nestrov_placer1.getMaxNesterovIter()
        || nestrov_iter2 >= nestrov_placer2.getMaxNesterovIter())
      break;
    // ONLY FOR DEBUGGING
    nestrov_placer1.updateDB();
    nestrov_placer2.updateDB();
    string file_name;
    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << i;
    ss >> file_name;;
    this->drawDies(file_name);
  }
  nestrov_placer1.updateDB();
  nestrov_placer2.updateDB();
}
}
