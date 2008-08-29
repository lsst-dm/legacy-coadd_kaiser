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
    unsigned int psfCols = psfKernel.getCols();
    unsigned int psfRows = psfKernel.getRows();
    unsigned int blurredPsfCols = 2 * psfCols - 1;
    unsigned int blurredPsfRows = 2 * psfRows - 1;
    // make image of psf, flip it in cols and rows, and zero-pad it
    // use imagePtr because Image.replaceSubImage requires a ptr
    lsst::afw::image::Image<double>::ImagePtrT reflPsfImagePtr(
        new lsst::afw::image::Image<double>(psfCols, psfRows));
    double dumImSum;
    psfKernel.computeImage(*reflPsfImagePtr, dumImSum, true);
    reflectImage(*reflPsfImagePtr);
    lsst::afw::image::Image<double> paddedReflPsfImage(blurredPsfCols, blurredPsfRows);
    vw::BBox2i centerBox(psfKernel.getCtrCol(), psfKernel.getCtrRow(), psfCols, psfRows);
    paddedReflPsfImage.replaceSubImage(centerBox, reflPsfImagePtr);
    
    throw std::runtime_error("afw does not yet support convolution with an Image");
    // convolve zero-padded reflected image of psf kernel with psf kernel
//    _blurredPsfImage = lsst::afw::math::convolve(paddedReflPsfImage, psfKernel, true);
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
    // this would be simpler with a proper STL iterator, but vw::PixelIterator is slow and read-only
    const unsigned int nCols = image.getCols();
    const unsigned int nRows = image.getRows();
    vw::MemoryStridingPixelAccessor<PixelType> frontAcc(image.getIVwPtr()->origin());
    vw::MemoryStridingPixelAccessor<PixelType> backAcc(image.getIVwPtr()->origin());
    backAcc.advance(nCols - 1, nRows - 1);
    unsigned int nColsToSwap = nCols;
    const unsigned int nRowsToSwap = (nRows + 1) / 2;
    const bool oddNRows = (nRows % 2 != 0);
    for (unsigned int row = 0; row < nRowsToSwap; ++row) {
        PixelType *frontPtr = &(*frontAcc);
        PixelType *backPtr = &(*backAcc);
        if (oddNRows && (row == nRowsToSwap - 1)) {
            nColsToSwap = nCols / 2;
        }
        for (unsigned int col = 0; col < nColsToSwap; ++col) {
            PixelType temp = *frontPtr;
            *frontPtr = *backPtr;
            *backPtr = temp;
        }
    }
}