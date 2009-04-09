// -*- LSST-C++ -*-
/**
* @brief Component of Kaiser coadd.
*
* @file
*
* @author Russell Owen
*/
#include <algorithm> // for swap

#include "lsst/pex/exceptions.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/coadd/kaiser.h"

namespace pexExcept = lsst::pex::exceptions;
namespace afwMath = lsst::afw::math;
namespace afwImage = lsst::afw::image;
namespace coaddKaiser = lsst::coadd::kaiser;

// local functions
namespace {

    /**
     * @brief Reflect an image in place.
     */
    void reflectImage(afwImage::Image<coaddKaiser::CoaddComponent::pixelType> &image) {
        typedef afwImage::Image<coaddKaiser::CoaddComponent::pixelType> ImageT;
        typedef ImageT::x_iterator XIterator;
        
        const unsigned int nRows = image.getHeight();
        const unsigned int nFullRowsToSwap = nRows / 2;
        const bool oddNRows = (nRows % 2 != 0);
         // use x_at(xLast, y) instead of row_end(y) to get the reverse row iterator
         // because row_end(y) starts one beyond the last pixel
        const int xLast = static_cast<int>(image.getWidth()) - 1;
        for (int yFwd = 0, yRev = image.getHeight() - 1; yFwd < nFullRowsToSwap; ++yFwd, --yRev) {
            for (XIterator fwdPtr = image.row_begin(yFwd), revPtr = image.x_at(xLast, yRev);
                fwdPtr != image.row_end(yFwd); ++fwdPtr, --revPtr) {
                std::swap(*fwdPtr, *revPtr);
            }
        }
        if (oddNRows) {
            const unsigned int yCtr = nFullRowsToSwap;
            const unsigned int halfCols = image.getWidth() / 2;
            XIterator const fwdEndCtr = image.x_at(halfCols + 1, yCtr); // +1 is for 1 beyond last pixel to swap
            XIterator fwdPtr = image.row_begin(yCtr);
            XIterator revPtr = image.x_at(xLast, yCtr);
            for ( ; fwdPtr != fwdEndCtr; ++fwdPtr, --revPtr) {
                std::swap(*fwdPtr, *revPtr);
            }
        }
    }

} // anonymous namespace

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
    _blurredExposure(scienceExposure.getMaskedImage().getWidth(), scienceExposure.getMaskedImage().getHeight()),
    _blurredPsfImage(psfKernel.getWidth() * 2 - 1, psfKernel.getHeight() * 2 - 1)
{
    computeSigmaSq(scienceExposure);
    computeBlurredPsf(psfKernel);
    computeBlurredExposure(scienceExposure, psfKernel);
};
        
void coaddKaiser::CoaddComponent::addToCoadd(ExposureD const &coadd) {
    throw LSST_EXCEPT(pexExcept::RuntimeErrorException, "Not implemented");
};

/**
 * @brief compute _sigmaSq
 *
 * @ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeSigmaSq(
    ExposureF const &scienceExposure    ///< science Exposure
) {
    typedef afwImage::MaskedImage<float, afwImage::MaskPixel>::x_iterator XIterator;
    typedef afwImage::MaskedImage<float, afwImage::MaskPixel>::y_iterator YIterator;

    MaskedImageF scienceMI(scienceExposure.getMaskedImage());
    // compute a vector containing only the good pixels, then take the median of that
    std::vector<double> varianceList(scienceMI.getHeight() * scienceMI.getWidth());
    std::vector<double>::iterator varIter = varianceList.begin();
    for (int y = 0; y != scienceMI.getHeight(); ++y) {
        for (XIterator ptr = scienceMI.row_begin(y); ptr != scienceMI.row_end(y); ++ptr) {
            if (ptr.mask() != 0) {
                continue;
            }
            *varIter = ptr.variance();
            ++varIter;
        }
    }
    // varIter now points to the end of the list
    _sigmaSq = coaddKaiser::medianBinapprox(varianceList.begin(), varIter);
// eventually something like the following will work directly on the variance image,
// but for now makeStatistics does not ignore masked pixels so is not usable; see PR #749
//     afwMath::Statistics varStats = afwMath::makeStatistics(*(scienceMI.getVariance()), afwMath::MEDIAN);
//     _sigmaSq = varStats.getValue(math::MEDIAN);
};
        
/**
 * @brief Compute _blurredPsfImage = psfKernel convolved with psfKernel(-r)
 *
 * @ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeBlurredPsf(
    afwMath::Kernel const &psfKernel    ///< PSF kernel
) {
    unsigned int psfWidth = psfKernel.getWidth();
    unsigned int psfHeight = psfKernel.getHeight();
    unsigned int blurredPsfWidth = 2 * psfWidth - 1;
    unsigned int blurredPsfHeight = 2 * psfHeight - 1;
    afwImage::Image<double> reflPsfImage(psfWidth, psfHeight);
    double dumImSum;
    psfKernel.computeImage(reflPsfImage, dumImSum, true);
    reflectImage(reflPsfImage);
    afwImage::Image<double> paddedReflPsfImage(blurredPsfWidth, blurredPsfHeight, 0);
    afwImage::BBox centerBox(afwImage::PointI(psfKernel.getCtrX(), psfKernel.getCtrY()), psfWidth, psfHeight);
    afwImage::Image<double>(paddedReflPsfImage, centerBox, false) <<= reflPsfImage;
    
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
    int edgeBit = scienceExposure.getMaskedImage().getMask()->getMaskPlane("EDGE");
    afwMath::convolve(_blurredExposure.getMaskedImage(), scienceExposure.getMaskedImage(), psfKernel, edgeBit, false);
    if (scienceExposure.hasWcs()) {
        _blurredExposure.setWcs(*scienceExposure.getWcs());
    }
};
