#ifndef INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#define INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#include "Parser.h"
#include "CImg.h"
namespace Color {
const unsigned char BLACK[] = {0, 0, 0};
const unsigned char DIM_GRAY[] = {105, 105, 105};
const unsigned char WHITE[] = {255, 255, 255};
const unsigned char BLUE[] = {0, 0, 255};
const unsigned char MINT[] = {3, 252, 219};
const unsigned char RED[] = {255, 0, 0};
const unsigned char PINK[] = {255, 51, 255};
const unsigned char LIGHT_YELLOW[] = {255, 236, 196};
const unsigned char LIGHT_GREEN[] = {73, 235, 52};
}

class Drawer {
  using Image = cimg_library::CImg<unsigned char>;
 public:
  Drawer(uint width, uint height);
  ~Drawer();
  void drawCell(int ll_x, int ll_y, int ur_x, int ur_y);
  void drawHybridBond(int ll_x, int ll_y, int ur_x, int ur_y);
  void saveImg(const string &file_name);
  void linkImg(Image *image);
  void setCellColor(const unsigned char *const color);
  Image *getImage() const;
  int getDieId() const;
  void setDieId(int die_id);

 private:
  int die_id_ = 0;
  uint width_ = 0;
  uint height_ = 0;
  string file_path = "../output/images/";
  Image *image_;
  bool is_linked = false;
  unsigned char color_[3] = {0, 0, 0};

};

#endif //INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
