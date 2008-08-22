#include <iostream>
#include <vector>
#include "lsst/coadd/kaiser.h"

int main() {
    typedef double dataType;
    const int nPts = 100;
    std::vector<dataType> dataArr(nPts);
    for (int i = 0; i < nPts; ++i) {
        dataArr[i] = static_cast<dataType>(i) / static_cast<dataType>(nPts);
    }
    double med = lsst::coadd::kaiser::medianBinapprox(dataArr.begin(), dataArr.end());
    std::cout << "Median = " << med << "; expected value is 0.5" << std::endl;
}