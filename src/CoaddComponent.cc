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
    _sigmaSq(0),
    _blurredExposure(),
    _blurredPsf()
{};

lsst::coadd::kaiser::CoaddComponent::CoaddComponent(
    Exposure const &scienceExposure,   ///< science Exposure
    lsst::afw::math::Kernel const &psfKernel    ///< PSF of science Exposure
) :
    _sigmaSq(0),
    _blurredExposure(),
    _blurredPsf()
{
    computeSigmaSq(scienceExposure);
//    computeBlurredPsf();
    computeBlurredExposure(scienceExposure, psfKernel);
};
        
void lsst::coadd::kaiser::CoaddComponent::addToCoadd(Exposure &coadd) {
    throw std::runtime_error("Not implemented");
};
        
void lsst::coadd::kaiser::CoaddComponent::computeSigmaSq(Exposure const &scienceExposure) {
    typedef lsst::afw::image::MaskedPixelAccessor<pixelType, lsst::afw::image::maskPixelType> ExposureAccessor;

    std::vector<pixelType> varianceList(scienceExposure.getCols() * scienceExposure.getRows());
    std::vector<pixelType>::iterator varIter = varianceList.begin();
    ExposureAccessor expRow(scienceExposure);
    for (unsigned int row = 0; row < scienceExposure.getRows(); ++row, expRow.nextRow()) {
        ExposureAccessor expCol = expRow;
        for (unsigned int col = 0; col < scienceExposure.getCols(); ++col, expCol.nextCol()) {
            if (*expCol.mask != 0) {
                continue;
            }
            *varIter = *expCol.variance;
            ++varIter;
        }
    }
    _sigmaSq = lsst::coadd::kaiser::medianBinapprox(varianceList.begin(), varIter);
};
        
/**
* @brief Compute _blurredPsf = psfKernel convolved with psfKernel(-r)
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
    int edgeBit = scienceExposure.getMask()->getMaskPlane("EDGE");
    _blurredExposure = lsst::afw::math::convolve(scienceExposure, psfKernel, edgeBit, false);
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