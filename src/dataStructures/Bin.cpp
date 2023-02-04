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
#include "Chip.h"

namespace VLSI_backend {
Chip::NestrovPlacer::Bin::Bin()
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

Chip::NestrovPlacer::Bin::Bin(int x, int y, int lx, int ly, int ux, int uy, float targetDensity)
    : Bin() {
  x_ = x;
  y_ = y;
  lx_ = lx;
  ly_ = ly;
  ux_ = ux;
  uy_ = uy;
  targetDensity_ = targetDensity;
}

Chip::NestrovPlacer::Bin::~Bin() {
  x_ = y_ = 0;
  lx_ = ly_ = ux_ = uy_ = 0;
  nonPlaceArea_ = instPlacedArea_ = fillerArea_ = nonPlaceAreaUnscaled_
      = instPlacedAreaUnscaled_ = 0;
  electroPhi_ = electroForceX_ = electroForceY_ = 0;
  density_ = targetDensity_ = 0;
}

const int64_t Chip::NestrovPlacer::Bin::binArea() const {
  return static_cast<int64_t>(dx()) * static_cast<int64_t>(dy());
}

float Chip::NestrovPlacer::Bin::density() const {
  return density_;
}

float Chip::NestrovPlacer::Bin::targetDensity() const {
  return targetDensity_;
}

float Chip::NestrovPlacer::Bin::electroForceX() const {
  return electroForceX_;
}

float Chip::NestrovPlacer::Bin::electroForceY() const {
  return electroForceY_;
}

float Chip::NestrovPlacer::Bin::electroPhi() const {
  return electroPhi_;
}

void Chip::NestrovPlacer::Bin::setDensity(float density) {
  density_ = density;
}

void Chip::NestrovPlacer::Bin::setTargetDensity(float density) {
  targetDensity_ = density;
}

void Chip::NestrovPlacer::Bin::setElectroForce(float electroForceX, float electroForceY) {
  electroForceX_ = electroForceX;
  electroForceY_ = electroForceY;
}

void Chip::NestrovPlacer::Bin::setElectroPhi(float phi) {
  electroPhi_ = phi;
}

int Chip::NestrovPlacer::Bin::x() const {
  return x_;
}

int Chip::NestrovPlacer::Bin::y() const {
  return y_;
}

int Chip::NestrovPlacer::Bin::lx() const {
  return lx_;
}

int Chip::NestrovPlacer::Bin::ly() const {
  return ly_;
}

int Chip::NestrovPlacer::Bin::ux() const {
  return ux_;
}

int Chip::NestrovPlacer::Bin::uy() const {
  return uy_;
}

int Chip::NestrovPlacer::Bin::cx() const {
  return (ux_ + lx_) / 2;
}

int Chip::NestrovPlacer::Bin::cy() const {
  return (uy_ + ly_) / 2;
}

int Chip::NestrovPlacer::Bin::dx() const {
  return (ux_ - lx_);
}

int Chip::NestrovPlacer::Bin::dy() const {
  return (uy_ - ly_);
}

void Chip::NestrovPlacer::Bin::setNonPlaceArea(int64_t area) {
  nonPlaceArea_ = area;
}

void Chip::NestrovPlacer::Bin::setNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ = area;
}

void Chip::NestrovPlacer::Bin::setInstPlacedArea(int64_t area) {
  instPlacedArea_ = area;
}

void Chip::NestrovPlacer::Bin::setInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ = area;
}

void Chip::NestrovPlacer::Bin::setFillerArea(int64_t area) {
  fillerArea_ = area;
}

void Chip::NestrovPlacer::Bin::addNonPlaceArea(int64_t area) {
  nonPlaceArea_ += area;
}

void Chip::NestrovPlacer::Bin::addInstPlacedArea(int64_t area) {
  instPlacedArea_ += area;
}

void Chip::NestrovPlacer::Bin::addNonPlaceAreaUnscaled(int64_t area) {
  nonPlaceAreaUnscaled_ += area;
}

void Chip::NestrovPlacer::Bin::addInstPlacedAreaUnscaled(int64_t area) {
  instPlacedAreaUnscaled_ += area;
}

void Chip::NestrovPlacer::Bin::addFillerArea(int64_t area) {
  fillerArea_ += area;
}

}