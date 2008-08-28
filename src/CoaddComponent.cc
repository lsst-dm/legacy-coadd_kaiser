// -*- LSST-C++ -*-
/**
* @brief Component of Kaiser coadd.
*
* @file
*
* @author Russell Owen
*/
#include "lsst/coadd/kaiser.h"

template<typename PixelType>
void reflectImage(lsst::afw::image::Image<PixelType> &image);

lsst::coadd::kaiser::CoaddComponent::CoaddComponent()
:
    LsstBase(typeid(this)),
    _sigmaSq(0),
    _blurredExposure(),
    _blurredPsfImage()
{};

lsst::coadd::kaiser::CoaddComponent::CoaddComponent(
    Exposure const &scienceExposure,   ///< science Exposure
    lsst::afw::math::Kernel const &psfKernel    ///< PSF of science Exposure
) :
    LsstBase(typeid(this)),
    _sigmaSq(0),
    _blurredExposure(),
    _blurredPsfImage()
{
    computeSigmaSq(scienceExposure);
//    computeBlurredPsf();
    computeBlurredExposure(scienceExposure, psfKernel);
};
        
void lsst::coadd::kaiser::CoaddComponent::addToCoadd(Exposure &coadd) {
    throw std::runtime_error("Not implemented");
};

        
void lsst::coadd::kaiser::CoaddComponent::computeSigmaSq(Exposure const &scienceExposure) {
    typedef lsst::afw::image::MaskedPixelAccessor<pixelType, lsst::afw::image::maskPixelType> MaskedPixelAccessor;

    MaskedImage scienceMI(scienceExposure.getMaskedImage());
    const unsigned int nCols(scienceMI.getCols());
    const unsigned int nRows(scienceMI.getRows());
    std::vector<pixelType> varianceList(nCols * nRows);
    std::vector<pixelType>::iterator varIter = varianceList.begin();
    MaskedPixelAccessor miRow(scienceMI);
    for (unsigned int row = 0; row < nRows; ++row, miRow.nextRow()) {
        MaskedPixelAccessor miCol = miRow;
        for (unsigned int col = 0; col < nCols; ++col, miCol.nextCol()) {
            if (*miCol.mask != 0) {
                continue;
            }
            *varIter = *miCol.variance;
            ++varIter;
        }
    }
    _sigmaSq = lsst::coadd::kaiser::medianBinapprox(varianceList.begin(), varIter);
};
        
/**
* @brief Compute _blurredPsfImage = psfKernel convolved with psfKernel(-r)
*/
void lsst::coadd::kaiser::CoaddComponent::computeBlurredPsf(
    lsst::afw::math::Kernel const &psfKernel    ///< PSF kernel
) {
    throw std::runtime_error("Not implemented");
};

/**
* @brief Compute _blurredEposure = scienceExposure convolved with psfKernel
*
* Warning: if you want scienceExposure convolved with psfKernel(-r)
* (the standard Kaiser thing to do) then feed in psfKernel(-r)
*/
void lsst::coadd::kaiser::CoaddComponent::computeBlurredExposure(
    Exposure const &scienceExposure,    ///< science exposure
    lsst::afw::math::Kernel const &psfKernel    ///< PSF kernel
) {
    // getMaskPlane should be a static function, but meanwhile...
    MaskedImage scienceMI(scienceExposure.getMaskedImage());
    int edgeBit = scienceMI.getMask()->getMaskPlane("EDGE");
    MaskedImage blurredMI(lsst::afw::math::convolve(scienceMI, psfKernel, edgeBit, false));
    _blurredExposure = Exposure(blurredMI, scienceExposure.getWcs());
};

/**
* @brief reflect an image in-place
*/
template<typename PixelType>
void reflectImage(lsst::afw::image::Image<PixelType> &image) {
    vw::PixelIterator<PixelType> front(*image.getVWPtr());
    vw::PixelIterator<PixelType> back(*image.getVWPtr(), image.getCols()-1, image.getRows()-1);
    unsigned long nToSwap = image.getCols() * image.getRows() / 2;
    for (unsigned long i = 0; i < nToSwap; ++i, ++front, --back) {
        PixelType temp = *front;
        *front = *back;
        *back = temp;
    }
}