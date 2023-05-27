#include "Drawer.h"

Drawer::Drawer(uint width, uint height) {
  width_ = width;
  height_ = height;
  image_ = new Image(width_, height_, 1, 3, 255);
}
Drawer::~Drawer() {
  if (is_linked)
    return;
  delete image_;
}
void Drawer::drawCell(int ll_x, int ll_y, int ur_x, int ur_y) {
  if (ll_x == ur_x)
    ur_x += 1;
  if (ll_y == ur_y)
    ur_y += 1;
  image_->draw_rectangle(ll_x, ll_y, ur_x, ur_y, Color::DIM_GRAY);
  image_->draw_rectangle(ll_x+1, ll_y+1, ur_x-1, ur_y-1, color_);
}
void Drawer::drawHybridBond(int ll_x, int ll_y, int ur_x, int ur_y) {
  int x = static_cast<int>((ll_x + ur_x) / 2);
  int y = static_cast<int>((ll_y + ur_y) / 2);
  image_->draw_circle(x, y, 1, Color::MINT);
}
void Drawer::saveImg(const string &file_name) {
  string save_file_name = file_path + file_name + ".bmp";
  image_->save_bmp(save_file_name.c_str());
}
void Drawer::linkImg(Drawer::Image *image) {
  delete image_;
  is_linked = true;
  image_ = image;
}
void Drawer::setCellColor(const unsigned char *const color) {
  color_[0] = color[0];
  color_[1] = color[1];
  color_[2] = color[2];
}
Drawer::Image *Drawer::getImage() const {
  return image_;
}
int Drawer::getDieId() const {
  return die_id_;
}
void Drawer::setDieId(int die_id) {
  die_id_ = die_id;
}
