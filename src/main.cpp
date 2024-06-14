#include "utils/globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"
#include "utils/Timer.hpp"

#ifdef HDF5_INTERFACE
#include <H5Cpp.h>
#endif

template<
    element_validity::uint n, 
    element_validity::uint s, 
    element_validity::uint p
>
void processData(
    element_validity::uint nNodesPerElem,
    element_validity::uint nElements,
    std::vector<element_validity::fp_t> nodes
) {
    element_validity::Timer timer;
    element_validity::ValidityChecker<n, s, p> checker;
    checker.setPrecisionTarget(.01);
    checker.setMaxSubdiv(30);
    const uint nCoordPerElem = nNodesPerElem*n*2;
    std::vector<element_validity::fp_t> results(nElements);
    std::vector<element_validity::fp_t> element(nCoordPerElem);

    for (uint e=1020; e<1100; ++e) {
        for(uint i=0; i<nCoordPerElem; ++i) {
            element.at(i) = nodes.at(e*nCoordPerElem + i);
        }
        timer.start();
        results.at(e) = checker.maxTimeStep(element);
        timer.stop();
        std::cout << "t = " << results.at(e) << "  ";
        std::cout << timer.read<std::chrono::microseconds>() << "us" << std::endl;
        timer.reset();
    }
}

int main(int argc, char** argv) {
    element_validity::Interval().init();

    if (argc <= 1) {
        std::cout << "Please provide filename" << std::endl;
        return 1;
    }

    uint dimension;
    uint nNodesPerElem;
    uint nElements;
    std::vector<element_validity::fp_t> nodes;

    {
        // Read data from HDF5 dataset
        std::cout << "Reading file..." << std::endl;
        std::string filename = argv[1];
        H5::H5File file(filename, H5F_ACC_RDONLY);

        const auto H5UINT = H5::PredType::NATIVE_INT;
        file.openDataSet("Dimension").read(&dimension, H5UINT);
        file.openDataSet("NumberOfHighOrderNodes").read(&nNodesPerElem, H5UINT);
        file.openDataSet("NumberOfSimplices").read(&nElements, H5UINT);

        H5::DataSet dataset = file.openDataSet("Nodes");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t nDataEntries = dataspace.getSimpleExtentNpoints();
        assert(4 * nElements * nNodesPerElem * dimension == nDataEntries);
        H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<char*> buffer(nDataEntries);
        dataset.read(buffer.data(), datatype, dataspace);
        nodes.reserve(nElements * nNodesPerElem * dimension * 2);
        for (hsize_t i=0; 2*i<nDataEntries; ++i) {
            element_validity::Rational rat;
            rat.num() = buffer[2*i];
            rat.den() = buffer[2*i+1];
            nodes.push_back(static_cast<element_validity::fp_t>(rat));
        }
        std::cout << "Done." << std::endl;
    }

    #define IFPROC(n, s, p, cp) if (dimension == n && nNodesPerElem == cp) \
        processData<n, s, p>(nNodesPerElem, nElements, nodes);
    IFPROC(1, 1, 1, 2)
    else IFPROC(1, 1, 2, 3)
    else IFPROC(1, 1, 3, 4)
    else IFPROC(1, 1, 4, 5)
    else IFPROC(1, 1, 5, 6)
    else IFPROC(2, 2, 1, 3)
    else IFPROC(2, 2, 2, 6)
    else IFPROC(2, 2, 3, 10)
    else IFPROC(3, 3, 1, 4)
    // IFPROC(3, 3, 2, 10)
    #undef IFPROC
    else throw std::invalid_argument("Not implemented");

    return 0;
}
