#include "db.h"
#include "lefin.h"
#include "defin.h"
#include "lefout.h"
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

void write_lef(odb::dbLib *lib, const std::string &path) {
  auto *logger = new utl::Logger(nullptr);
  std::ofstream os;
  os.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  os.open(path.c_str());
  odb::lefout writer(logger, os);
  writer.writeTechAndLib(lib);
}
void makeShirkedLef(odb::dbDatabase *db_database_original, const std::string &new_lef_path) {
  auto libs_original = db_database_original->getLibs();
  auto tech_original = db_database_original->getTech();

  auto db_database = odb::dbDatabase::create();
  auto tech = odb::dbTech::create(db_database);
  tech->setLefUnits(tech_original->getLefUnits());
  tech->setManufacturingGrid(tech_original->getManufacturingGrid());
  tech->setLefVersion(tech_original->getLefVersion());
  tech->setDbUnitsPerMicron(tech_original->getDbUnitsPerMicron());

  for (auto layer_original : tech_original->getLayers()) {
    odb::dbTechLayer::create(tech, layer_original->getName().c_str(), layer_original->getType());
  }

  int lib_number = 0;

  for (auto lib_original : libs_original) {
    lib_number++;
    auto lib = odb::dbLib::create(db_database, lib_original->getConstName());
    lib->setLefUnits(lib_original->getLefUnits());
    for (auto site_original : lib_original->getSites()) {
      auto site = odb::dbSite::create(lib, site_original->getConstName());
      site->setWidth(site_original->getWidth());
      site->setHeight(site_original->getHeight());
    }

    for (auto master_original : lib_original->getMasters()) {
      int master_origin_x, master_origin_y;
      odb::dbMaster *master = odb::dbMaster::create(lib, master_original->getConstName());
      master->setType(master_original->getType());
      if (master_original->getEEQ())
        master->setEEQ(master_original->getEEQ());
      if (master_original->getLEQ())
        master->setLEQ(master_original->getLEQ());
      master->setWidth(master_original->getWidth());
      master->setHeight(master_original->getHeight());

      master_original->getOrigin(master_origin_x, master_origin_y);
      master->setOrigin(master_origin_x, master_origin_y);

      auto site = lib->findSite(master_original->getSite()->getConstName());
      master->setSite(site);

      if (master_original->getSymmetryX())
        master->setSymmetryX();
      if (master_original->getSymmetryY())
        master->setSymmetryY();
      if (master_original->getSymmetryR90())
        master->setSymmetryR90();

      master->setMark(master_original->isMarked());
      master->setSequential(master_original->isSequential());
      master->setSpecialPower(master_original->isSpecialPower());

      for (auto m_term_original : master_original->getMTerms()) {
        auto db_m_term = odb::dbMTerm::create(master, m_term_original->getConstName(),
                                              m_term_original->getIoType(),
                                              m_term_original->getSigType(),
                                              m_term_original->getShape());

        for (auto pin_original : m_term_original->getMPins()) {
          auto db_m_pin = odb::dbMPin::create(db_m_term);
          for (auto geometry : pin_original->getGeometry()) {
            odb::dbBox::create(db_m_pin,
                               geometry->getTechLayer(),
                               geometry->xMin(), geometry->yMin(),
                               geometry->xMax(), geometry->yMax());
          }
        }
      }
      master->setFrozen();
    }
    write_lef(lib, new_lef_path + std::to_string(lib_number));
    write_lef(lib_original, "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input_original.lef");
  }

}
int main() {
  std::string lef_path = "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input.lef";
  std::string new_lef_path = "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input_new.lef";
  std::string def_path = "../test/benchmarks/standard/ispd/ispd18_test1/ispd18_test1.input.def";
  std::string db_path = "../output/dbFiles/test.db";

  odb::dbDatabase *db_database = odb::dbDatabase::create();
  parseLef(db_database, lef_path);
  makeShirkedLef(db_database, new_lef_path);

}

