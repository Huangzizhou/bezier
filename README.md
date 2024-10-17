# High-Order Continuous Geometrical Validity
## Installation steps
The code requires a x86 compliant processor, a Linux machine and the gcc compiler.
We will extend to other operating systems and architectures after publication.

```
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

HDF5 datasets from the benchmark are in `meshes.zip`.