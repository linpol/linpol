#include <cmath>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "io-xyz.h"
#include "xyz-water-utils.h"

#include "mbpol.h"

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: test-mbpol water.xyz"
                  << std::endl;
        return 0;
    }

    std::cout << std::scientific << std::setprecision(9);

    std::vector<std::string> elements;
    std::vector<double> crd;

    try {
        std::ifstream ifs(argv[1]);

        if (!ifs)
            throw std::runtime_error("could not open the XYZ file");

        std::string comment;
        kit::io::load_xyz(ifs, comment, elements, crd);
    } catch (const std::exception& e) {
        std::cerr << " ** Error ** : " << e.what() << std::endl;
        return 1;
    }

    if (!kit::is_water(elements, crd, false)) {
        std::cerr << " ** Error ** : not water?" << std::endl;
        return 1;
    }

    const size_t nw = elements.size()/3;

    x2o::mbpol pot;

    double grd[9*nw];
    double E = pot(nw, &(crd[0]), grd);

    std::cout << "       E = " << E << "\n"
              << "E[nogrd] = " << pot(nw, &(crd[0])) << "\n"
              << std::endl;

    const double eps = 1.0e-4;

    for (int n = 0; n < 9*nw; ++n) {
        double tmp[9*nw];
        const double x_orig = crd[n];

        crd[n] = x_orig + eps;
        const double Ep = pot(nw, &(crd[0]), tmp);

        crd[n] = x_orig - eps;
        const double Em = pot(nw, &(crd[0]), tmp);

        const double gfd = 0.5*(Ep - Em)/eps;
        crd[n] = x_orig;

//        std::cout << grd[n] << ' ' << gfd << '\n';
        std::cout << n << ' ' << std::fabs(grd[n] - gfd) << '\n';
    }

    return 0;
}
