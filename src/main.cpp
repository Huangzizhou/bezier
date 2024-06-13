#include "utils/globals.hpp"
#include "numbers/Interval.hpp"
#include "utils/combinatorics.hpp"
#include "validity/element_validity.hpp"
#include "utils/Timer.hpp"

#ifdef HDF5_INTERFACE
#include <H5Cpp.h>
#endif

int main(int argc, char** argv) {
    element_validity::Interval().init();
    using element_validity::fp_t;
    using element_validity::Validity;
    using element_validity::ValidityChecker;
    // switch () {
    // case Validity::valid: std::cout << "valid"; break;
    // case Validity::invalid: std::cout << "invalid"; break;
    // case Validity::uncertain: std::cout << "uncertain"; break;
    // }

    uint dimension;
    uint nNodesPerElem;
    uint nElements;
    std::vector<fp_t> nodes;

    if (argc > 1) {
        std::cout << "Reading file..." << std::endl;
        // Read data from HDF5 dataset
        std::string filename = argv[1];
        H5::H5File file(filename, H5F_ACC_RDONLY);

        const auto H5UINT = H5::PredType::NATIVE_INT;
        file.openDataSet("Dimension").read(&dimension, H5UINT);
        file.openDataSet("NumberOfHighOrderNodes").read(&nNodesPerElem, H5UINT);
        file.openDataSet("NumberOfSimplices").read(&nElements, H5UINT);

		std::vector<std::string> components;
        H5::DataSet dataset = file.openDataSet("Nodes");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t nDataEntries = dataspace.getSimpleExtentNpoints();
        assert(4 * nElements * nNodesPerElem * dimension == nDataEntries);
        H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<char*> buffer(nDataEntries);
        dataset.read(buffer.data(), datatype, dataspace);
        for (hsize_t i = 0; i < nDataEntries; ++i)
            components.emplace_back(buffer[i]);

        nodes.reserve(nElements * nNodesPerElem * dimension * 2);


        for (uint i=0; 2*i<nDataEntries; ++i) {
            element_validity::Rational rat;
            rat.num() = components.at((2*i)+0);
            rat.den() = components.at((2*i)+1);
            nodes.emplace_back(rat);
		}
        std::cout << "Done." << std::endl;
    }


    using element_validity::Timer;
    Timer timer;
    ValidityChecker<2, 2, 2> v222;
    v222.setPrecisionTarget(.01);
    v222.setMaxSubdiv(30);
    const uint nCoordPerElem = nNodesPerElem*dimension*2;
    std::vector<fp_t> results(nElements);
    std::vector<fp_t> element(nCoordPerElem);

    // timer.start();
    // const double t = v222.maxTimeStep({
    //     0,0, 0,0,
    //     1,1, 0,0,
    //     .15,.85, .7,.7,
    //     .5,.5, -.1,-.1,
    //     .57,.57, .46,.25,
    //     .43,.43, .25,.46
    // });
    // timer.stop();
    // std::cout << "test t = " << t << std::endl;
    // std::cout << timer.read<std::chrono::microseconds>() << "us" << std::endl;
    // timer.reset();

    // for (uint e=0; e<nElements; ++e) {
    for (uint e=1020; e<1090; ++e) {
        for(uint i=0; i<nCoordPerElem; ++i) {
            element.at(i) = nodes.at(e*nCoordPerElem + i);
        }
        timer.start();
        results.at(e) = v222.maxTimeStep(element);
        timer.stop();
        std::cout << "t = " << results.at(e) << "  ";
        std::cout << timer.read<std::chrono::microseconds>() << "us" << std::endl;
        timer.reset();
    }
    
    return 0;
}
