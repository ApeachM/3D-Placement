#include "Drawer.h"

Drawer::Drawer(uint width, uint height, uint margin_x, uint margin_y) {
  width_ = width;
  height_ = height;
  margin_x_ = margin_x;
  margin_y_ = margin_y;
  image_ = Image(width_ + 2 * margin_x_, height_ + margin_y_ * 2, 1, 3, 255);
  image_.draw_rectangle(margin_x_, margin_y_, width_ + margin_x_, height_ + margin_y_, COLOR::BLACK);
  image_.draw_rectangle(margin_x_ + 1, margin_y_ + 1, width_ + margin_x_ - 1, height_ + margin_y_ - 1, COLOR::WHITE);
}
Drawer::~Drawer() {
}
void Drawer::drawCell(int ll_x, int ll_y, int ur_x, int ur_y, int DIE_ID) {
  if (ll_x == ur_x) {
    ur_x += 2;
    ll_x -= 1;
  }
  if (ll_y == ur_y) {
    ur_y += 2;
    ll_y -= 1;
  }
  ll_x += static_cast<int>(margin_x_ - average_cell_width_ / 2);
  ll_y += static_cast<int>(margin_y_ - average_cell_height_ / 2);
  ur_x += static_cast<int>(margin_x_ - average_cell_width_ / 2);
  ur_y += static_cast<int>(margin_y_ - average_cell_height_ / 2);
  image_.draw_rectangle(ll_x, ll_y, ur_x, ur_y, COLOR::DIM_GRAY);
  if (DIE_ID == 0)
    image_.draw_rectangle(ll_x + 1, ll_y + 1, ur_x - 1, ur_y - 1, default_cell_color_);
  else if (DIE_ID == 1)
    image_.draw_rectangle(ll_x + 1, ll_y + 1, ur_x - 1, ur_y - 1, top_cell_color_);
  else if (DIE_ID == 2)
    image_.draw_rectangle(ll_x + 1, ll_y + 1, ur_x - 1, ur_y - 1, bottom_cell_color_);
}
void Drawer::drawHybridBond(int ll_x, int ll_y, int ur_x, int ur_y) {
  int x = static_cast<int>((ll_x + ur_x) / 2);
  int y = static_cast<int>((ll_y + ur_y) / 2);
  x += static_cast<int>(margin_x_ - average_cell_width_ / 2);
  y += static_cast<int>(margin_y_ - average_cell_height_ / 2);
  image_.draw_circle(x, y, 2, COLOR::BLACK);
  image_.draw_circle(x, y, 1, COLOR::RED);
}
void Drawer::saveImg(const string &file_name) {
  string save_file_name = file_path + file_name + ".bmp";
  image_.save_bmp(save_file_name.c_str());
}
void Drawer::setTopCellColor(const unsigned char *const color) {
  top_cell_color_[0] = color[0];
  top_cell_color_[1] = color[1];
  top_cell_color_[2] = color[2];
}
void Drawer::setMarginX(uint margin_x) {
  margin_x_ = margin_x;
}
void Drawer::setMarginY(uint margin_y) {
  margin_y_ = margin_y;
}
void Drawer::setAverageCellWidth(uint average_cell_width) {
  average_cell_width_ = average_cell_width;
}
void Drawer::setAverageCellHeight(uint average_cell_height) {
  average_cell_height_ = average_cell_height;
}
