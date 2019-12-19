# eigen-osqp

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![CI](https://github.com/jrl-umi3218/eigen-osqp/workflows/CI%20of%20eigen-osqp/badge.svg?branch=master)](https://github.com/jrl-umi3218/eigen-osqp/actions?query=workflow%3A%22CI+of+eigen-osqp%22)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](http://jrl-umi3218.github.io/eigen-osqp/doxygen/HEAD/index.html)

eigen-osqp allows to use the OSQP solver with the Eigen3 library.

It supports dense and sparse matrix.

## Installation

### Ubuntu LTS (16.04, 18.04, 20.04)

```bash
# Make sure you have required tools
sudo apt install apt-transport-https lsb-release ca-certificates gnupg
# Add our key
sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key 892EA6EE273707C6495A6FB6220D644C64666806
# Add our repository (stable versions)
sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-release $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Use this to setup the HEAD version
# sudo sh -c 'echo "deb https://dl.bintray.com/gergondet/multi-contact-head $(lsb_release -sc) main" | sudo tee /etc/apt/sources.list.d/multi-contact.list'
# Update packages list
sudo apt update
# Install packages
sudo apt install libeigen-osqp-dev
```

### From source

#### Dependencies

To compile you need the following tools:

* C++ compiler with C++11 support
* [CMake]() >= 3.1
* [Doxygen](http://www.stack.nl/~dimitri/doxygen/): to generate documentation (optional)
* [Boost](http://www.boost.org) >= 1.49
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
* [OSQP](https://github.com/oxfordcontrol/osqp) >= 0.6.0

#### Building

```sh
mkdir _build
cd _build
cmake [options] ..
make && make intall
```
