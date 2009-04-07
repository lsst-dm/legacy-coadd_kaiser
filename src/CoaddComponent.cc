// -*- LSST-C++ -*-
/**
* @brief Component of Kaiser coadd.
*
* @file
*
* @author Russell Owen
*/
#include <algorithms> // for swap

#include "lsst/pex/exceptions.h"
#include "lsst/coadd/kaiser.h"

template<typename PixelType>

namespace pexExcept = lsst::pex::exceptions;
namespace afwMath = lsst::afw::math;
namespace afwImage = lsst::afw::image;
namespace coaddKaiser = coaddKaiser;

/**
 * @brief CoaddComponent constructor
 *
 * @todo: handle asymmetric kernels properly. scienceExposure should be convolved with * the *reflected* PSF,
 * but this requires significant extra work (perhaps new convolution functions or kernels) to handle
 * spatially varying kernels. For now, for expediency, I allow spatially varying kernels but convolve with
 * the un-reflected PSF.
 * @todo: address Robert Lupton's concerns about handling the background. He feels we should not
 * subtract the background from science exposures before adding them to the template, but I fail to see
 * how we can avoid doing so. Otherwise the background of the template will vary pixel by pixel
 * depending on how many good pixels from the various science exposures contributed to a give pixel
 * of the template.
 *
 * @ingroup coadd::kaiser
 */ 
coaddKaiser::CoaddComponent::CoaddComponent(
    ExposureF const &scienceExposure,   ///< background-subtracted science Exposure
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
        
void coaddKaiser::CoaddComponent::addToCoadd(ExposureD &coadd) {
    throw pexExcept::Runtime("Not implemented");
};

/**
 * @brief compute _sigmaSq
 *
 * @ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeSigmaSq(
    ExposureF const &scienceExposure    ///< science Exposure
) {
    typedef afwImage::MaskedPixelAccessor<float, afwImage::MaskPixel> MaskedPixelAccessorF;
    typedef typename afwImage::MaskedImage<float, afwImage::MaskPixel>::x_iterator x_iteratorF;
    typedef typename afwImage::MaskedImage<float, afwImage::MaskPixel>::y_iterator y_iteratorF;

    MaskedImageF scienceMI(scienceExposure.getMaskedImage());
    std::vector<double> varianceList(scienceMI.getHeight() * scienceMI.getWidth());
    std::vector<double>::iterator varIter = varianceList.begin();
    for (int y = 0; y != scienceMI.getHeight(); ++y) {
        for (x_iteratorF ptr = scienceExposure.row_begin(y); ptr != scienceExposure.row_end(y); ++ptr) {
            if (ptr->mask() != 0) {
                continue;
            }
            *varIter = ptr->variance();
            ++varIter;

        }
    }
    _sigmaSq = coaddKaiser::medianBinapprox(varianceList.begin(), varIter);
};
        
/**
 * @brief Compute _blurredPsfImage = psfKernel convolved with psfKernel(-r)
 *
 * @ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeBlurredPsf(
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
    afwImage::BBox centerBox(afwImage::PointI(psfKernel.getCtrCol(), psfKernel.getCtrRow()), psfCols, psfRows);
    afwImage::Image<double>(paddedReflPsfImage, centerBox, false) << reflPsfImage;
    
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
void coaddKaiser::CoaddComponent::computeBlurredExposure(
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
void coaddKaiser::reflectImage(afwImage::Image<PixelType> &image) {
    typedef afw::Image::Image<PixelType> ImageT;
    
    const unsigned int nFullRowsToSwap = nRows / 2;
    const bool oddNRows = (nRows % 2 != 0);
     // use x_at(xLast, y) instead of row_end(y) to get the reverse row iterator
     // because row_end(y) starts one beyond the last pixel
    const int xLast = image.getWidth() - 1;
    for (int yFwd = 0, yRev = image.getHeight() - 1; yFwd < nFullRowsToSwap; ++yFwd, --yRev) {
        for (ImageT::x_iterator fwdPtr = image.row_begin(yFwd), revPtr = image.x_at(xLast, yRev);
            fwdPtr != image.row_end(yFwd); ++fwdPtr, --revPtr) {
            std::swap(*fwdPtr, *revPtr);
        }
    }
    if (oddNRows) {
        const unsigned int yCtr = nFullRowsToSwap;
        const unsigned int halfCols = image.getWidth() / 2;
        ImageT::x_iterator fwdEndCtr = image.x_at(halfCols + 1, yCtr); // +1 for 1 beyond last pixel to swap
        for (ImageT::x_iterator fwdPtr = image.row_begin(yCtr), revPtr = image.x_at(xLast, yCtr);
            fwdPtr != fwdEndCtr; ++fwdPtr, --revPtr) {
            std::swap(*fwdPtr, *revPtr);
        }
    }
}