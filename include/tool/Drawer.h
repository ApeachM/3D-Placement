#ifndef INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#define INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#include "Parser.h"
#include "CImg.h"
namespace COLOR {
const unsigned char BLACK[] = {0, 0, 0};
const unsigned char DIM_GRAY[] = {105, 105, 105};
const unsigned char WHITE[] = {255, 255, 255};
const unsigned char BLUE[] = {0, 0, 255};
const unsigned char MINT[] = {3, 252, 219};
const unsigned char RED[] = {255, 0, 0};
const unsigned char PINK[] = {255, 51, 255};
const unsigned char LIGHT_YELLOW[] = {255, 236, 196};
const unsigned char LIGHT_GREEN[] = {73, 235, 52};
const unsigned char LIGHT_RED[] = {255, 150, 150};
const unsigned char LIGHT_MINT[] = {70, 150, 150};
}

class Drawer {
  using Image = cimg_library::CImg<unsigned char>;
 public:
  Drawer(uint width, uint height, uint margin_x, uint margin_y);
  ~Drawer();
  void drawCell(int ll_x, int ll_y, int ur_x, int ur_y, int DIE_ID);
  void drawHybridBond(int ll_x, int ll_y, int ur_x, int ur_y);
  void saveImg(const string &file_name);
  void setTopCellColor(const unsigned char *const color);
  void setBottomCellColor(const unsigned char *const color) {
    bottom_cell_color_[0] = color[0];
    bottom_cell_color_[1] = color[1];
    bottom_cell_color_[2] = color[2];
  }
  void setMarginX(uint margin_x);
  void setMarginY(uint margin_y);
  void setAverageCellWidth(uint average_cell_width);
  void setAverageCellHeight(uint average_cell_height);

 private:
  int die_id_ = 0;
  uint width_ = 0;
  uint height_ = 0;
  uint margin_x_ = 0;
  uint margin_y_ = 0;
  uint average_cell_width_ = 0;
  uint average_cell_height_ = 0;
  string file_path = "../output/images/";
  Image image_;
  unsigned char default_cell_color_[3] = {0, 0, 0};
  unsigned char top_cell_color_[3] = {0, 0, 0};
  unsigned char bottom_cell_color_[3] = {0, 0, 0};

};

#endif //INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
