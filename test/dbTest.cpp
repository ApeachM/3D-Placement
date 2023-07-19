#include "db.h"
#include "lefin.h"
#include "defin.h"
#include <iostream>

void parseLef(odb::dbDatabase *db_database, const std::string &lef_path) {
  utl::Logger logger_;
  odb::lefin lef_reader(db_database, &logger_, false);
  odb::dbLib *lib = lef_reader.createTechAndLib("lib_name", lef_path.c_str());
  odb::dbTech *tech = db_database->getTech();

  if (lib != nullptr || tech != nullptr) {
    std::cout << "Lef parsing is succeed." << std::endl;
  } else {
    std::cout << "Lef parsing is failed." << std::endl;
  }
}
void parseDef(odb::dbDatabase *db_database, const std::string &def_path) {
  utl::Logger logger_;
  odb::defin def_reader(db_database, &logger_);
  std::vector<odb::dbLib *> search_libs;
  for (odb::dbLib *lib : db_database->getLibs())
    search_libs.push_back(lib);
  odb::dbChip *chip = def_reader.createChip(search_libs, def_path.c_str());
  if (chip) {
    std::cout << "Def parsing is succeed." << std::endl;
  } else {
    std::cout << "Def parsing is failed." << std::endl;
  }
}
void saveDb(odb::dbDatabase *db_database, const std::string &db_path) {
  FILE *stream = std::fopen(db_path.c_str(), "w");
  if (stream) {
    db_database->write(stream);
    std::fclose(stream);
  }
}

void dbExample(odb::dbDatabase *db_database) {

  std::cout << db_database->getChip()->getBlock()->getBBox()->getDX() << std::endl;

  odb::dbBlock *block = db_database->getChip()->getBlock();
  for (int i = 0; i < 4; ++i) {
    std::cout << std::endl;
  }
  // instance cloning into circuit variables
  odb::dbSet<odb::dbInst> db_instances = block->getInsts();
  int cellIdx = 0;

  for (odb::dbInst *db_instance : db_instances) {
    if (cellIdx % 10 == 0) {
      std::cout << "db_instance->getMaster()->getLib()->getDbUnitsPerMicron(): "
                << db_instance->getMaster()->getLib()->getDbUnitsPerMicron() << std::endl;
      int x, y;
      auto type = db_instance->getMaster()->getType();
      if (!type.isCore() && !type.isBlock()) continue;

      std::cout << "db_instance->getName(): " << db_instance->getName() << std::endl;
      db_instance->getOrigin(x, y);
      std::cout << "getOrigin: " << x << " " << y << std::endl;
      std::cout << "db_instance->getMaster()->getName(): " << db_instance->getMaster()->getName() << std::endl;
      std::cout << "db_instance->getMaster()->getHeight(): " << db_instance->getMaster()->getHeight() << std::endl;
      std::cout << "db_instance->getMaster()->getWidth(): " << db_instance->getMaster()->getWidth() << std::endl;
      std::cout << "db_instance->getMaster()->getLib()->getName(): " << db_instance->getMaster()->getLib()->getName()
                << std::endl << std::endl;

      for (auto masterTerminal : db_instance->getMaster()->getMTerms()) {
        std::cout << std::endl;
        std::cout << "terminal name: " << masterTerminal->getName() << std::endl;
        std::cout << "masterTerminal->getSigType().getString(): " << masterTerminal->getSigType().getString()
                  << std::endl;
        std::cout << "masterTerminal->getIoType().getString(): " << masterTerminal->getIoType().getString()
                  << std::endl;
        std::cout << "masterTerminal->getIndex(): " << masterTerminal->getIndex() << std::endl;
        std::cout << "for (auto pin: masterTerminal->getMPins())" << std::endl;
        for (auto pin : masterTerminal->getMPins()) {
          for (auto box : pin->getGeometry()) {
            std::cout << "\t" << box->xMin() << " " << box->xMax() << "  " << box->yMin() << " " << box->yMax()
                      << std::endl;
          }
        }
      }
      std::cout << std::endl;

      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        if (instanceTerminal->getNet())
          std::cout << instanceTerminal->getNet()->getName() << std::endl;
        else
          std::cout << "None net." << std::endl;
        std::cout << " " << instanceTerminal->getMTerm()->getName()
                  << " " << instanceTerminal->getIoType().getString() << " " << x << " " << y << std::endl;;
      }

      std::cout << std::endl << std::endl;
      std::cout << "get origin: " << x << " " << y << std::endl;
      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        std::cout << x << " " << y << std::endl;
      }
      db_instance->setOrigin(100, 200);
      db_instance->getOrigin(x, y);
      std::cout << "get origin: " << x << " " << y << std::endl;
      for (auto instanceTerminal : db_instance->getITerms()) {
        instanceTerminal->getAvgXY(&x, &y);
        std::cout << x << " " << y << std::endl;
      }

      std::cout << "instance end." << std::endl << std::endl;
      cellIdx++;
    }
  }
/*

  std::cout << std::endl << std::endl << std::endl << std::endl;
  std::cout << "OpenDB Tutorial starts." << std::endl;

  dbBlock *block = parser_.pseudo_db_database_->getChip()->getBlock();
  dbSet<dbInst> db_instances = block->getInsts();

  // How to access all instances in database, and print the instance name
  int numOfInstance = 0;
  for (dbInst *db_instance : db_instances) {
    if (!db_instance->getName().empty())
      numOfInstance++;
  }
  std::cout << "The number of Instances: " << numOfInstance << std::endl;

  // Let's access only the first one instance, and explore some other methods.
  dbInst *db_inst = *db_instances.begin();
  std::cout << "The name of the instance: " << db_inst->getName() << std::endl;
  // As the library of the instance, we can know the library name, and cell width and height.
  dbMaster *lib = db_inst->getMaster();
  std::cout << "The name of the library of the instance: " << lib->getName() << std::endl;
  std::cout << "The width of the instance: " << db_inst->getMaster()->getWidth() << std::endl;
  std::cout << "The height of the instance: " << db_inst->getMaster()->getHeight() << std::endl << std::endl;

  std::cout << "The pin libraries of the instance" << std::endl;
  for (auto library_terminals : db_inst->getMaster()->getMTerms()) {
    std::cout << " ";
    std::cout << "Name: " << library_terminals->getName() << "\t";
    std::cout << "Signal Type: " << library_terminals->getSigType().getString() << "\t";
    std::cout << "Terminal Index: " << library_terminals->getIndex() << "\t";
    std::cout << std::endl;
  }
  std::cout << "The number of the terminals of instance: " << db_inst->getMaster()->getMTermCount() << std::endl;
  std::cout << std::endl;

  std::cout << "The pin data" << std::endl;
  for (auto instance_terminal : db_inst->getITerms()) {
    std::cout << " ";
    std::cout << "Terminal Name: " << instance_terminal->getMTerm()->getName() << "\t";
    std::cout << "Signal Type: " << instance_terminal->getSigType().getString() << "\t";
    std::cout << "Terminal Index: " << instance_terminal->getMTerm()->getIndex() << std::endl;
    std::cout << "  Connected Net: ";
    if (instance_terminal->getNet() != nullptr) {
      std::cout << instance_terminal->getNet()->getName() << std::endl;
      if (instance_terminal->getNet()->getName() != "clk") {
        std::cout << "   connected Pins correspond to the net" << std::endl;
        for (auto terminals_in_net : instance_terminal->getNet()->getITerms()) {
          std::cout << "   "
               << terminals_in_net->getMTerm()->getName() << " in "
               << terminals_in_net->getInst()->getName() << "("
               << terminals_in_net->getInst()->getMaster()->getName() << ") instance" << std::endl;
        }
      }
    } else {
      std::cout << "None" << std::endl;
    }
  }

  // # 2. traverse nets
  for (int i = 0; i < 3; ++i) {
    std::cout << std::endl;
  }
  // access the nets in the circuit
  std::cout << "Signal Type\tIoType\tTerminal Name" << std::endl;
  for (auto net : block->getNets()) {
    for (auto terminal : net->getITerms()) {
    }
    for (auto terminal : net->getBTerms()) {
      std::cout << terminal->getSigType().getString() << "\t";
      std::cout << terminal->getIoType().getString() << "\t";
      std::cout << terminal->getName() << std::endl;
      for (auto bPin : terminal->getBPins()) {
        std::cout << " bPin ID: " << bPin->getId() << std::endl;
        for (auto box : bPin->getBoxes()) {
          std::cout << "  box xMin: " << box->xMin() << std::endl;
          std::cout << "  box xMax: " << box->xMax() << std::endl;
          std::cout << "  box DX: " << box->getDX() << std::endl;
          std::cout << "  box yMin: " << box->yMin() << std::endl;
          std::cout << "  box yMax: " << box->yMax() << std::endl;
          std::cout << "  box DY: " << box->getDY() << std::endl;
          std::cout << "  box width: " << box->getWidth() << std::endl;
          std::cout << "  box dir: " << box->getDir() << std::endl;
          if (box->getBlockVia())*/
/* meaningless in this case*//*

            std::cout << "  via name: " << box->getBlockVia()->getName() << std::endl;
        }

      }
    }
  }


  // # 3. Get Die information
  odb::Rect Die, Core;
  block->getDieArea(Die);
  block->getCoreArea(Core);
  std::cout << Die.dx() << " " << Die.dy() << std::endl;
  std::cout << Core.dx() << " " << Core.dy() << std::endl;
*/

}
int main() {
  std::string lef_path = "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input.lef";
  std::string def_path = "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input.def";
  std::string db_path = "../output/dbFiles/test.db";

}

