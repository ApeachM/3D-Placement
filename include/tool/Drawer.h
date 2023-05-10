#ifndef INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#define INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
#include "Parser.h"
#include "CImg.h"
namespace Color {
const unsigned char BLACK[] = {0, 0, 0};
const unsigned char DIM_GRAY[] = {105, 105, 105};
const unsigned char WHITE[] = {255, 255, 255};
const unsigned char BLUE[] = {3, 252, 219};
const unsigned char RED[] = {255, 0, 0};
const unsigned char PINK[] = {255, 51, 255};
const unsigned char LIGHT_YELLOW[] = {255, 236, 196};
const unsigned char LIGHT_GREEN[] = {73, 235, 52};
}

class Drawer {
  using Image = cimg_library::CImg<unsigned char>;
 public:
  Drawer(uint width, uint height) {
    width_ = width;
    height_ = height;
    image_ = Image(width_, height_, 1, 3, 255);
  }
  void drawRect(int ll_x, int ll_y, int ur_x, int ur_y, const unsigned char *color) {
    image_.draw_rectangle(ll_x, ll_y, ur_x, ur_y, color);
  }
  void saveImg(string file_name) {
    file_name = file_path + file_name;
    image_.save_bmp(file_name.c_str());
  }

  int getDieId() const {
    return die_id_;
  }
  void setDieId(int die_id) {
    die_id_ = die_id;
  }

 private:
  int die_id_ = 0;
  uint width_ = 0;
  uint height_ = 0;
  string file_path = "../output/images/";
  Image image_;

};

#endif //INC_3D_PLACEMENT_INCLUDE_TOOL_DRAWER_H_
