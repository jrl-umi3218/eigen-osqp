# eigen-osqp

[![License](https://img.shields.io/badge/License-BSD%202--Clause-green.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![CI](https://github.com/jrl-umi3218/eigen-osqp/workflows/CI%20of%20eigen-osqp/badge.svg?branch=master)](https://github.com/jrl-umi3218/eigen-osqp/actions?query=workflow%3A%22CI+of+eigen-osqp%22)
[![Documentation](https://img.shields.io/badge/doxygen-online-brightgreen?logo=read-the-docs&style=flat)](http://jrl-umi3218.github.io/eigen-osqp/doxygen/HEAD/index.html)

eigen-osqp allows to use the OSQP solver with the Eigen3 library.

It supports dense and sparse matrix.

## Installation

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
