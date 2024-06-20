#include "utils/globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"
#include "utils/Timer.hpp"
#include "utils/Settings.hpp"

#include <memory>
#include <fstream>

#ifdef HDF5_INTERFACE
#include <H5Cpp.h>

constexpr char SEP = ',';

template<
    element_validity::uint n, 
    element_validity::uint s, 
    element_validity::uint p
>
void processData(
    element_validity::Settings args,
    element_validity::uint nNodesPerElem,
    element_validity::uint nElements,
    std::vector<element_validity::fp_t> nodes,
    std::ostream *out = nullptr
) {
    element_validity::Timer timer;
    element_validity::ValidityChecker<n, s, p> checker;
    checker.setPrecisionTarget(.01);
    // checker.setMaxSubdiv(30);
    const uint nCoordPerElem = nNodesPerElem*n*2;
    std::vector<element_validity::fp_t> results;
    const uint ne = std::min(nElements, args.lastElem - args.firstElem);
    results.reserve(ne);

    if (out)
        *out << "ID" << SEP
            << "max_time_step" << SEP
            << "space_depth" << SEP
            << "time_depth" << SEP
            << "iterations" << SEP
            << "microseconds" << std::endl;
    std::vector<element_validity::fp_t> element(nCoordPerElem);
    for (uint e=0; e<ne; ++e) {
        for(uint i=0; i<nCoordPerElem; ++i) {
            element.at(i) = nodes.at((e + args.firstElem)*nCoordPerElem + i);
        }
	    std::array<uint, 3> info;
        timer.start();
        const element_validity::fp_t t = checker.maxTimeStep(element, &info);
        timer.stop();
        results.push_back(t);
        if (out)
            *out << e + args.firstElem << SEP
                << element_validity::fp_fmt << t << SEP
                << std::get<0>(info) << SEP
                << std::get<1>(info) << SEP
                << std::get<2>(info) << SEP
                << timer.read<std::chrono::microseconds>() << std::endl;
        timer.reset();
    }
}

int main(int argc, char** argv) {
    element_validity::Interval().init();
    element_validity::Settings args(argc, argv);

    uint dimension;
    uint nNodesPerElem;
    uint nElements;
    std::vector<element_validity::fp_t> nodes;

    {
        // Read rational data from HDF5 dataset
        std::cout << "Reading file..." << std::endl;
        H5::H5File file(args.filePath, H5F_ACC_RDONLY);

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

    std::unique_ptr<std::ofstream> out;
    if (args.resultsPath.size() > 0)
        out = std::make_unique<std::ofstream>(args.resultsPath);
    #define IFPROC(n, s, p, cp) if (dimension == n && nNodesPerElem == cp) \
        processData<n, s, p>(args, nNodesPerElem, nElements, nodes, out.get());
    IFPROC(1, 1, 1, 2)
    else IFPROC(1, 1, 2, 3)
    else IFPROC(1, 1, 3, 4)
    else IFPROC(1, 1, 4, 5)
    else IFPROC(1, 1, 5, 6)
    else IFPROC(2, 2, 1, 3)
    else IFPROC(2, 2, 2, 6)
    // else IFPROC(2, 2, 3, 10)
    // else IFPROC(2, 2, 4, 15)
    else IFPROC(3, 3, 1, 4)
    // IFPROC(3, 3, 2, 10)
    #undef IFPROC
    else throw std::invalid_argument("Not implemented");

    return 0;
}
#else
#warning "HDF5 interface disabled, ignoring main"
int main(int argc, char** argv) { return 0; }
#endif