# EDA API

[![Build Status](https://img.shields.io/badge/develop-ongoing%20-green)]()

This repository is API for EDA development.
This API is based on [OpenDB](https://github.com/The-OpenROAD-Project/OpenDB) API.

Currently, this project is building for giving template of course work term project of POSTECH in Korea:

The purpose of term project is **Implementing placer**.

## External Dependencies

You can use `Dockerfile` in the `OpenDB` submodule,
or you are required to install `cmake`, `siwg`, `spdlog`, `boost` in the `Ubuntu` Environment.

## How to build

```shell
git clone --recurse-submodules https://github.com/ApeachM/EDA-API.git
cd EDA-API
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