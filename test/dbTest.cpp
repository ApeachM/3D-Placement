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
odb::dbDatabase *loadDB(const std::string &db_path) {
  odb::dbDatabase *db_database = odb::dbDatabase::create();
  std::ifstream file;
  file.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ios::eofbit);
  file.open(db_path, std::ios::binary);
  db_database->read(file);
  return db_database;
}
void write_lef(odb::dbLib *lib, const std::string &path) {
  auto *logger = new utl::Logger(nullptr);
  std::ofstream os;
  os.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  os.open(path.c_str());
  odb::lefout writer(logger, os);
  writer.writeTechAndLib(lib);
}
void makeShrankLef(odb::dbDatabase *db_database_original, const std::string &new_lef_path,
                   float shirk_ratio = 0.7, const std::string &which_die = "") {
  auto libs_original = db_database_original->getLibs();
  auto tech_original = db_database_original->getTech();

  auto db_database = odb::dbDatabase::create();
  auto tech = odb::dbTech::create(db_database);
  tech->setLefUnits(tech_original->getLefUnits());
  tech->setManufacturingGrid(tech_original->getManufacturingGrid());
  tech->setLefVersion(tech_original->getLefVersion());
  tech->setDbUnitsPerMicron(tech_original->getDbUnitsPerMicron());

  for (auto layer_original : tech_original->getLayers()) {
    odb::dbTechLayer::create(tech, (layer_original->getName() + which_die).c_str(), layer_original->getType());
  }

  int lib_number = 0;

  for (auto lib_original : libs_original) {
    lib_number++;
    auto lib = odb::dbLib::create(db_database, (lib_original->getName() + which_die).c_str());
    lib->setLefUnits(lib_original->getLefUnits());
    for (auto site_original : lib_original->getSites()) {
      auto site = odb::dbSite::create(lib, (site_original->getName() + which_die).c_str());
      // apply shrink factor on site
      site->setWidth(static_cast<int>(static_cast<float>(site_original->getWidth()) * shirk_ratio));
      site->setHeight(static_cast<int> (static_cast<float> (site_original->getHeight()) * shirk_ratio));
    }

    for (auto master_original : lib_original->getMasters()) {
      odb::dbMaster *master = odb::dbMaster::create(lib, (master_original->getName() + which_die).c_str());
      master->setType(master_original->getType());
      if (master_original->getEEQ())
        master->setEEQ(master_original->getEEQ());
      if (master_original->getLEQ())
        master->setLEQ(master_original->getLEQ());
      int width = static_cast<int>(static_cast<float>(master_original->getWidth()) * shirk_ratio);
      int height = static_cast<int>(static_cast<float>(master_original->getHeight()) * shirk_ratio);
      master->setWidth(width);
      master->setHeight(height);

      int master_origin_x, master_origin_y;
      master_original->getOrigin(master_origin_x, master_origin_y);
      master_origin_x = static_cast<int>(static_cast<float>(master_origin_x) * shirk_ratio);
      master_origin_y = static_cast<int>(static_cast<float>(master_origin_y) * shirk_ratio);
      master->setOrigin(master_origin_x, master_origin_y);

      std::string site_name = master_original->getSite()->getName() + which_die;
      auto site = lib->findSite(site_name.c_str());
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
            int x1, y1, x2, y2;
            x1 = static_cast<int>(static_cast<float>(geometry->xMin()) * shirk_ratio);
            y1 = static_cast<int>(static_cast<float>(geometry->yMin()) * shirk_ratio);
            x2 = static_cast<int>(static_cast<float>(geometry->xMax()) * shirk_ratio);
            y2 = static_cast<int>(static_cast<float>(geometry->yMax()) * shirk_ratio);
            odb::dbBox::create(db_m_pin, geometry->getTechLayer(), x1, y1, x2, y2);
          }
        }
      }
      master->setFrozen();
    }
    assert(lib_number == 1);
    write_lef(lib, new_lef_path);
  }
}
void makeTopDieLef(odb::dbDatabase *db_database, const std::string &lef_path) {
  // makeShrankLef(db_database, lef_path, 0.7, "_top");
  makeShrankLef(db_database, lef_path, 0.7);
}
void makeBottomDieLef(odb::dbDatabase *db_database, const std::string &lef_path) {
  // makeShrankLef(db_database, lef_path, 1.0, "_bottom");
  makeShrankLef(db_database, lef_path, 1.0);
}
void testRatioPartition(){
  auto db_database = loadDB("/home/mik077/tmp/tmp.yGF2jgykCL2/submodules/OpenROAD/src/par/test/gcd_0x100.db");

  int num_1 = 0;
  int num_2 = 0;
  for (auto inst : db_database->getChip()->getBlock()->getInsts()) {
    if (odb::dbIntProperty::find(inst, "partition_id")->getValue() == 0) {
      num_1++;
    } else if (odb::dbIntProperty::find(inst, "partition_id")->getValue() == 1) {
      num_2++;
    }
  }

  // total db inst number
  std::cout << "total db inst number: " << db_database->getChip()->getBlock()->getInsts().size() << std::endl;

  std::cout << "num_1: " << num_1 << std::endl;
  std::cout << "num_2: " << num_2 << std::endl;

  // partition ratio
  std::cout << "partition ratio (1): " << static_cast<float>(num_1) / static_cast<float>(num_1 + num_2) << std::endl;
  std::cout << "partition ratio (2): " << static_cast<float>(num_2) / static_cast<float>(num_1 + num_2) << std::endl;
}
void makeShrunkLef(){
  std::string dir_path = "../test/benchmarks/standard/ispd/ispd18_test1/";
  std::string design_name = "ispd18_test1";
  std::string lef_name = "ispd18_test1.input.lef";
  std::string def_name = "ispd18_test1.input.def";
  std::string db_path = "../output/dbFiles/test.db";

  odb::dbDatabase *db_database = odb::dbDatabase::create();
  parseLef(db_database, dir_path + design_name + ".input.lef");
  makeTopDieLef(db_database, dir_path + design_name + ".input_top.lef");
  makeBottomDieLef(db_database, dir_path + design_name + ".input_bottom.lef");
}
int main() {
}

