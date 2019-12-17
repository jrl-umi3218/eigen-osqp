# eigen-osqp

eigen-osqp allows to use the osqp QP solver with the Eigen3 library.

## Installing

### Manual

#### Dependencies

To compile you need the following tools:

 * [CMake]() >= 2.8
 * [pkg-config]()
 * [g++]()
 * [Boost](http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html) >= 1.49
 * [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) >= 3.2
 * [OSQP](https://github.com/oxfordcontrol/osqp) >= 0.6.0

#### Building

```sh
mkdir _build
cd _build
cmake [options] ..
make && make intall
```

Where the main options are:

 * `-DCMAKE_BUIlD_TYPE=Release` Build in Release mode
 * `-DCMAKE_INSTALL_PREFIX=some/path/to/install` default is `/usr/local`

