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
#include "NesterovPlacer.h"

namespace VLSI_backend {
Chip::NesterovPlacer::Bin::Bin()
    : x_(0),
      y_(0),
      lx_(0),
      ly_(0),
      ux_(0),
      uy_(0),
      nonPlaceArea_(0),
      instPlacedArea_(0),
      instPlacedAreaUnscaled_(0),
      nonPlaceAreaUnscaled_(0),
      fillerArea_(0),
      density_(0),
      targetDensity_(0),
      electroPhi_(0),
      electroForceX_(0),
      electroForceY_(0) {
}

Chip::NesterovPlacer::Bin::Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity)
    : Bin() {
  x_ = x;
  y_ = y;
  lx_ = lx;
  ly_ = ly;
  ux_ = ux;
  uy_ = uy;
  targetDensity_ = targetDensity;
}

Chip::NesterovPlacer::Bin::~Bin() {
  x_ = y_ = 0;
  lx_ = ly_ = ux_ = uy_ = 0;
  nonPlaceArea_ = instPlacedArea_ = fillerArea_ = nonPlaceAreaUnscaled_
      = instPlacedAreaUnscaled_ = 0;
  electroPhi_ = electroForceX_ = electroForceY_ = 0;
  density_ = targetDensity_ = 0;
}

const int64_t Chip::NesterovPlacer::Bin::binArea() const {
  return static_cast<int64_t>(dx()) * static_cast<int64_t>(dy());
}

float Chip::NesterovPlacer::Bin::density() const {
  return density_;
}

float Chip::NesterovPlacer::Bin::targetDensity() const {
  return targetDensity_;
}

float Chip::NesterovPlacer::Bin::electroForceX() const {
  return electroForceX_;
}

float Chip::NesterovPlacer::Bin::electroForceY() const {
  return electroForceY_;
}

float Chip::NesterovPlacer::Bin::electroPhi() const {
  return electroPhi_;
}

void Chip::NesterovPlacer::Bin::setDensity(float density) {
  density_ = density;
}

void Chip::NesterovPlacer::Bin::setTargetDensity(float density) {
  targetDensity_ = density;
}

void Chip::NesterovPlacer::Bin::setElectroForce(float electroForceX, float electroForceY) {
  electroForceX_ = electroForceX;
  electroForceY_ = electroForceY;
}

void Chip::NesterovPlacer::Bin::setElectroPhi(float phi) {
  electroPhi_ = phi;
}

int Chip::NesterovPlacer::Bin::x() const {
  return x_;
}

int Chip::NesterovPlacer::Bin::y() const {
  return y_;
}

int Chip::NesterovPlacer::Bin::lx() const {
  return lx_;
}

int Chip::NesterovPlacer::Bin::ly() const {
  return ly_;
}

int Chip::NesterovPlacer::Bin::ux() const {
  return ux_;
}

int Chip::NesterovPlacer::Bin::uy() const {
  return uy_;
}

int Chip::NesterovPlacer::Bin::cx() const {
  return (ux_ + lx_) / 2;
}

int Chip::NesterovPlacer::Bin::cy() const {
  return (uy_ + ly_) / 2;
}

int Chip::NesterovPlacer::Bin::dx() const {
  return (ux_ - lx_);
}

int Chip::NesterovPlacer::Bin::dy() const {
  return (uy_ - ly_);
}

void Chip::NesterovPlacer::Bin::setNonPlaceArea(int64_t area) {
  nonPlaceArea_ = area;
}

void Chip::NesterovPlacer::Bin::setNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ = area;
}

void Chip::NesterovPlacer::Bin::setInstPlacedArea(int64_t area) {
  instPlacedArea_ = area;
}

void Chip::NesterovPlacer::Bin::setInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ = area;
}

void Chip::NesterovPlacer::Bin::setFillerArea(int64_t area) {
  fillerArea_ = area;
}

void Chip::NesterovPlacer::Bin::addNonPlaceArea(int64_t area) {
  nonPlaceArea_ += area;
}

void Chip::NesterovPlacer::Bin::addInstPlacedArea(int64_t area) {
  instPlacedArea_ += area;
}

void Chip::NesterovPlacer::Bin::addNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ += area;
}

void Chip::NesterovPlacer::Bin::addInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ += area;
}

void Chip::NesterovPlacer::Bin::addFillerArea(int64_t area) {
  fillerArea_ += area;
}

}