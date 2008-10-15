// -*- LSST-C++ -*-
/**
* @brief Component of Kaiser coadd.
*
* @file
*
* @author Russell Owen
*/
#include "lsst/pex/exceptions.h"
#include "lsst/coadd/kaiser.h"

template<typename PixelType>
void reflectImage(lsst::afw::image::Image<PixelType> &image);

namespace pexExcept = lsst::pex::exceptions;
namespace afwMath = lsst::afw::math;
namespace afwImage = lsst::afw::image;

/**
 * @brief CoaddComponent constructor
 *
 * @todo: decide if we want to support spatially varying kernels. If we want both that
 * and scienceExposure to be convolved with the *reflected* PSF then we some significant work to do;
 * it will probably require new convolution functions. For now, for expediency, I allow spatially
 * varying kernel but do not convolve with the reflected PSF.
 *
 * @ingroup coadd::kaiser
 */ 
lsst::coadd::kaiser::CoaddComponent::CoaddComponent(
    ExposureF const &scienceExposure,   ///< science Exposure
    afwMath::Kernel const &psfKernel    ///< PSF of science Exposure
) :
    LsstBase(typeid(this)),
    _sigmaSq(0),
    _blurredExposure(scienceExposure.getMaskedImage().getCols(), scienceExposure.getMaskedImage().getRows()),
    _blurredPsfImage(psfKernel.getCols() * 2 - 1, psfKernel.getRows() * 2 - 1)
{
    computeSigmaSq(scienceExposure);
//    computeBlurredPsf();
    computeBlurredExposure(scienceExposure, psfKernel);
};
        
void lsst::coadd::kaiser::CoaddComponent::addToCoadd(ExposureD &coadd) {
    throw pexExcept::Runtime("Not implemented");
};

/**
 * @brief compute _sigmaSq
 *
 * @ingroup coadd::kaiser
 */
void lsst::coadd::kaiser::CoaddComponent::computeSigmaSq(
    ExposureF const &scienceExposure    ///< science Exposure
) {
    typedef afwImage::MaskedPixelAccessor<float, afwImage::maskPixelType> MaskedPixelAccessorF;

    MaskedImageF scienceMI(scienceExposure.getMaskedImage());
    const unsigned int nCols(scienceMI.getCols());
    const unsigned int nRows(scienceMI.getRows());
    std::vector<double> varianceList(nCols * nRows);
    std::vector<double>::iterator varIter = varianceList.begin();
    MaskedPixelAccessorF miRow(scienceMI);
    for (unsigned int row = 0; row < nRows; ++row, miRow.nextRow()) {
        MaskedPixelAccessorF miCol = miRow;
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
 *
 * @ingroup coadd::kaiser
 */
void lsst::coadd::kaiser::CoaddComponent::computeBlurredPsf(
    afwMath::Kernel const &psfKernel    ///< PSF kernel
) {
    unsigned int psfCols = psfKernel.getCols();
    unsigned int psfRows = psfKernel.getRows();
    unsigned int blurredPsfCols = 2 * psfCols - 1;
    unsigned int blurredPsfRows = 2 * psfRows - 1;
    afwImage::Image<double> reflPsfImage(psfCols, psfRows);
    double dumImSum;
    psfKernel.computeImage(reflPsfImage, dumImSum, true);
    reflectImage(reflPsfImage);
    afwImage::Image<double> paddedReflPsfImage(blurredPsfCols, blurredPsfRows);
    vw::BBox2i centerBox(psfKernel.getCtrCol(), psfKernel.getCtrRow(), psfCols, psfRows);
    paddedReflPsfImage.replaceSubImage(centerBox, reflPsfImage);
    
    // convolve zero-padded reflected image of psf kernel with psf kernel
    afwMath::convolve(_blurredPsfImage, paddedReflPsfImage, psfKernel, true);
};

/**
 * @brief Compute _blurredEposure = scienceExposure convolved with psfKernel
 *
 * @warning If you want scienceExposure convolved with psfKernel(-r)
 * (the standard Kaiser thing to do) then feed in psfKernel(-r)
 *
 * @ingroup coadd::kaiser
 */
void lsst::coadd::kaiser::CoaddComponent::computeBlurredExposure(
    ExposureF const &scienceExposure,   ///< science exposure
    afwMath::Kernel const &psfKernel    ///< PSF kernel
) {
    // getMaskPlane should be a static function, but meanwhile...
    MaskedImageF scienceMI(scienceExposure.getMaskedImage());
    int edgeBit = scienceMI.getMask()->getMaskPlane("EDGE");
    MaskedImageD blurredMI = _blurredExposure.getMaskedImage();
    afwMath::convolve(blurredMI, scienceMI, psfKernel, edgeBit, false);
    _blurredExposure.setWcs(scienceExposure.getWcs());
};

/**
 * @brief reflect an image in-place
 *
 * @ingroup coadd::kaiser
 */
template<typename PixelType>
void reflectImage(afwImage::Image<PixelType> &image) {
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