# 3D Placement with D2D Vertical Connections

[![Build Status](https://img.shields.io/badge/develop-ongoing%20-green)]()

This repository is 3D placer in VLSI back-end design, based on [ICCAD2022 Contest Problem B](http://iccad-contest.org).

This is based on [OpenDB](https://github.com/The-OpenROAD-Project/OpenDB) API.



## External Dependencies

You can use `Dockerfile` in the `OpenDB` submodule, or you are required to install `cmake`, `siwg`, `spdlog`, `boost` in the `Ubuntu` Environment.

You can use `etc/DependencyInstaller.sh` in linux system.

## How to build

```shell
git clone --recurse-submodules https://github.com/ApeachM/3D-Placement.git
cd 3D-Placement
git submodule update --force --recursive --init --remote
```

```shell
mkdir build
cd build
cmake ..
make
./placer
```



## Contributor

[Minjae Kim](https://github.com/ApeachM) in POSTECH
Email: kmj0824@postech.ac.kr
