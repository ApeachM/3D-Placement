# 3D Placement with D2D Vertical Connections

[![Build Status](https://img.shields.io/badge/develop-ongoing%20-green)]()

This is working on progress to be integrated into OpenROAD Project.

This repository is 3D placer in VLSI back-end design, based on [ICCAD2022 Contest Problem B](http://iccad-contest.org).

This is based on [OpenDB](https://github.com/The-OpenROAD-Project/OpenDB) API.



## External Dependencies

You can use [`Dockerfile`](submodules/OpenROAD/etc/DockerHelper.sh) in the [`OpenROAD`](submodules/OpenROAD) submodule, 
or you are required to install `cmake`, `siwg`, `spdlog`, `boost` in the `Ubuntu` Environment.

The shell script in [`etc/DependencyInstaller.sh`](etc/DependencyInstaller.sh) can give you help.

## How to build

```shell
git clone --recurse-submodules https://github.com/ApeachM/3D-Placement.git
cd 3D-Placement
```

```shell
mkdir build
cd build
cmake ..
make
./placer3D
```



## Contributor

[Minjae Kim](https://github.com/ApeachM) in POSTECH
Email: kmj0824@postech.ac.kr

[SoonHyun Kwon](https://github.com/kwonsh01) in POSTECH
Email: shyun010302@postech.ac.kr
