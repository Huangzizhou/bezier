#include "utils/globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"
#include "utils/combinatorics.hpp"
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
    using namespace element_validity;
    Timer timer;
    ValidityChecker<n, s, p> checker;
    checker.setPrecisionTarget(args.precision);
    checker.setMaxSubdiv(args.maxIterations);
    const uint nCoordPerElem = nNodesPerElem*n*2;
    std::vector<fp_t> results;
    const uint ne = args.numElem == 0 ? nElements : args.numElem;
    results.reserve(ne);

    if (out)
        *out << "ID" << SEP
            << "max_time_step" << SEP
            << "time_of_inversion" << SEP
            << "space_depth" << SEP
            << "microseconds" << SEP
            << "hierarchy" << SEP
            << "description" << std::endl;
    std::vector<fp_t> element(nCoordPerElem);
    for (uint e=0; e<ne; ++e) {
        for(uint i=0; i<nCoordPerElem; ++i) {
            element.at(i) = nodes.at((e + args.firstElem)*nCoordPerElem + i);
        }
        std::vector<uint> h;
        CheckerInfo info;
        fp_t tInv;
        timer.start();
        const fp_t t = checker.maxTimeStep(element, &h, 1, &tInv, &info);
        timer.stop();
        results.push_back(t);
        if (out) {
            *out << e + args.firstElem << SEP;
            *out << fp_fmt << t << SEP;
            *out << fp_fmt << tInv << SEP;
            *out << info.spaceDepth << SEP;
            *out << timer.read<std::chrono::microseconds>() << SEP;
            for (uint u : h) *out << u << ' ';
            *out << SEP;
            *out << info.description() << std::endl;
        }
        timer.reset();
    }
}

int main(int argc, char** argv) {
    using namespace element_validity;
    Interval().init();
    Settings args(argc, argv);

    uint dimension;
    uint nNodesPerElem;
    uint nElements;
    std::vector<fp_t> nodes;

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
            Rational rat;
            rat.num() = buffer[2*i];
            rat.den() = buffer[2*i+1];
            nodes.push_back(static_cast<fp_t>(rat));
        }
        std::cout << "Done." << std::endl;
    }

    std::unique_ptr<std::ostream> out;
    if (args.resultsPath.size() > 0)
        out = std::make_unique<std::ofstream>(args.resultsPath);
    std::ostream *const outptr = out ? out.get() : &std::cout;
    #define IFPROC(n, s, p) \
        if (dimension == n && nNodesPerElem == nControlGeoMap(n,s,p)) \
        processData<n, s, p>(args, nNodesPerElem, nElements, nodes, outptr);
    IFPROC(1, 1, 1)
    else IFPROC(1, 1, 2)
    else IFPROC(1, 1, 3)
    else IFPROC(1, 1, 4)
    else IFPROC(1, 1, 5)
    else IFPROC(2, 1, 1)
    else IFPROC(2, 1, 2)
    else IFPROC(2, 2, 1)
    else IFPROC(2, 2, 2)
    else IFPROC(2, 2, 3)
    else IFPROC(2, 2, 4)
    else IFPROC(3, 3, 1)
    else IFPROC(3, 3, 2)
    else IFPROC(3, 3, 3)
    else IFPROC(3, 3, 4)
    #undef IFPROC
    else throw std::invalid_argument("Not implemented");

    return 0;
}
#else
#warning "HDF5 interface disabled, ignoring main"
int main(int argc, char** argv) { return 0; }
#endif