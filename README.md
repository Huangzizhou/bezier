# High-Order Continuous Geometrical Validity
## Installation steps
```
git clone https://gitlab.com/fsichetti/robust-bezier-subdivision.git
cd robust-bezier-subdivision
mkdir build
cd build
cmake ..
make -j
```

Run `./tests/unit_tests` to check if the library works.
Run `./tests/eigen_tests` to check if the Eigen interface works.
Basic usage is `./bin yourdata.hdf5 -o youroutput.csv`.
Run the binary with `-?` to see instructions on all available flags.