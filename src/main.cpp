#include "element_validity.hpp"
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
element_validity::fp_t processData(
    element_validity::Settings args,
    element_validity::uint nNodesPerElem,
    element_validity::uint nElements,
    std::vector<element_validity::fp_t> nodes,
    std::ostream *out = nullptr
) {
    using namespace element_validity;
    Timer timer;
    ContinuousValidator<n, s, p> checker(args.numThreads);
    checker.setPrecisionTarget(args.precision);
    checker.setMaxSubdiv(args.maxIterations);
    StaticValidator<n, s, p> sChecker(args.numThreads);
    sChecker.setMaxSubdiv(args.preCheckMaxIter);
    const uint nCoordPerElem = nNodesPerElem*n*2;
    const uint lastElem =
        args.numElem == 0 ? nElements : args.firstElem + args.numElem;

    if (out)
        *out << "ID" << SEP
            << "max_time_step" << SEP
            << "time_of_inversion" << SEP
            << "space_depth" << SEP
            << "time_depth" << SEP
            << "microseconds" << SEP
            << "hierarchy" << SEP
            << "description" << std::endl;

    // Continuous check
    fp_t minT = 1;
    for (uint e=args.firstElem; e<lastElem; ++e) {
        const uint elemOffset = e*nCoordPerElem;
        span<fp_t> element(nodes.data() + elemOffset, nCoordPerElem);
        std::vector<uint> h;
        Validator::Info info;
        fp_t tInv;
        if (args.preCheck) {
            const Validity v0 = sChecker.isValidStart(element);
            if (v0 != Validity::valid) {
                if (v0 == Validity::uncertain)
                    std::cerr <<
                        "Warning: static check could not determine whether \
                        element " << e << " is valid at t=0" << std::endl;
                continue;
            }
        }
        timer.start();
        const fp_t t = checker.maxTimeStep(element, &h, nullptr, &tInv, &info);
        timer.stop();
        minT = std::min(minT, t);
        const double microseconds =
            static_cast<double>(timer.read<std::chrono::nanoseconds>()) / 1000;
        if (out) {
            *out << e << SEP;
            *out << fp_fmt << t << SEP;
            *out << fp_fmt << tInv << SEP;
            *out << info.spaceDepth << SEP;
            *out << info.timeDepth << SEP;
            *out << microseconds << SEP;
            for (uint u : h) *out << u << ' ';
            *out << SEP;
            *out << info.description() << std::endl;
        }
        timer.reset();
    }
    return minT;
}

int main(int argc, char** argv) {
    using namespace element_validity;
    Interval::init();
    Settings args(argc, argv);
    if (args.abort) return 0;

    uint dimension;
    uint nNodesPerElem;
    uint nElements;
    std::vector<fp_t> nodes;

    {
        // Read rational data from HDF5 dataset
        std::cout << "Reading file..." << std::flush;
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
        std::cout << " Done." << std::endl;
    }

    std::cout << "Processing..." << std::flush;
    std::unique_ptr<std::ostream> out;
    if (args.resultsPath.size() > 0)
        out = std::make_unique<std::ofstream>(args.resultsPath);
    std::ostream *const outptr = out ? out.get() : &std::cout;
    fp_t mtt = 0;
    #define IFPROC(n, s, p) \
        if (dimension == n && nNodesPerElem == nControlGeoMap(n,s,p)) \
        mtt = processData<n, s, p>(\
            args, nNodesPerElem, nElements, nodes, outptr);
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
    else IFPROC(3, 1, 1)
    else IFPROC(3, 2, 1)
    else IFPROC(3, 3, 1)
    else IFPROC(3, 3, 2)
    else IFPROC(3, 3, 3)
    // else IFPROC(3, 3, 4)
    #undef IFPROC
    else throw std::invalid_argument("Not implemented");
    std::cout << " Done." << std::endl;
    std::cout << "Max time step: " << mtt << std::endl;

    return 0;
}
#else
#warning "HDF5 interface disabled, ignoring main"
int main(int argc, char** argv) { return 0; }
#endif
