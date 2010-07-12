// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
* \brief Component of Kaiser coadd.
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
     * \brief Reflect an image in place.
     */
    void reflectImage(afwImage::Image<coaddKaiser::CoaddComponent::pixelType> &image) {
        typedef coaddKaiser::CoaddComponent::ImageCC::x_iterator XIteratorCC;
        
        if ((image.getHeight() < 1) || (image.getWidth() < 1)) {
            throw LSST_EXCEPT(pexExcept::RangeErrorException, "Image width and/or height is 0");
        }
        
        const int nRows = static_cast<int>(image.getHeight());
        const int nFullRowsToSwap = nRows / 2;
        const bool isOddNRows = (nRows % 2 != 0);
         // use x_at(xLast, y) instead of row_end(y) to get the reverse row iterator
         // because row_end(y) starts one beyond the last pixel
        const int xLast = static_cast<int>(image.getWidth()) - 1;
        for (int yFwd = 0, yRev = nRows - 1; yFwd < nFullRowsToSwap; ++yFwd, --yRev) {
            for (XIteratorCC fwdPtr = image.row_begin(yFwd), revPtr = image.x_at(xLast, yRev);
                fwdPtr != image.row_end(yFwd); ++fwdPtr, --revPtr) {
                std::swap(*fwdPtr, *revPtr);
            }
        }
        if (isOddNRows) {
            const unsigned int yCtr = nFullRowsToSwap;
            const unsigned int halfCols = image.getWidth() / 2;
            XIteratorCC const fwdEndCtr = image.x_at(halfCols + 1, yCtr); // +1 for 1 after last pixel to swap
            XIteratorCC fwdPtr = image.row_begin(yCtr);
            XIteratorCC revPtr = image.x_at(xLast, yCtr);
            for ( ; fwdPtr != fwdEndCtr; ++fwdPtr, --revPtr) {
                std::swap(*fwdPtr, *revPtr);
            }
        }
    }

} // anonymous namespace

/**
 * \brief CoaddComponent constructor
 *
 * \todo: handle asymmetric kernels properly. scienceExposure should be convolved with the *reflected* PSF,
 * but this requires significant extra work (perhaps new convolution functions or kernels) to handle
 * spatially varying kernels. For now, for expediency, I allow spatially varying kernels but convolve with
 * the un-reflected PSF.
 *
 * \ingroup coadd::kaiser
 */ 
coaddKaiser::CoaddComponent::CoaddComponent(
    ExposureF const &scienceExposure,   ///< science Exposure with the background subtracted
    lsst::afw::math::Kernel const &psfKernel,   ///< PSF of science Exposure
    bool normalizePsf                   ///< normalize psf
) :
    lsst::daf::base::Citizen(typeid(this)),
    _sigmaSq(0),
    _blurredExposure(scienceExposure.getWidth(), scienceExposure.getHeight()),
    _blurredPsfImage(psfKernel.getWidth() * 2 - 1, psfKernel.getHeight() * 2 - 1, 0),
    _normalizePsf(normalizePsf)
{
    computeSigmaSq(scienceExposure);
    computeBlurredPsf(psfKernel);
    computeBlurredExposure(scienceExposure, psfKernel);
};

/**
 * \brief compute _sigmaSq
 *
 * \ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeSigmaSq(
    ExposureF const &scienceExposure    ///< science Exposure
) {
    typedef ExposureF::MaskedImageT::x_iterator XIteratorF;

    ExposureF::MaskedImageT scienceMI = scienceExposure.getMaskedImage();
    // compute a vector containing only the good pixels, then take the median of that
    std::vector<double> varianceList(scienceMI.getHeight() * scienceMI.getWidth());
    std::vector<double>::iterator varIter = varianceList.begin();
    for (int y = 0, yEnd = scienceMI.getHeight(); y < yEnd; ++y) {
        for (XIteratorF ptr = scienceMI.row_begin(y); ptr != scienceMI.row_end(y); ++ptr) {
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
 * \brief Compute _blurredPsfImage = psfKernel convolved with psfKernel(-r)
 *
 * \ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeBlurredPsf(
    lsst::afw::math::Kernel const &psfKernel    ///< PSF kernel
) {
    int const psfWidth = psfKernel.getWidth();
    int const psfHeight = psfKernel.getHeight();
    int const paddedWidth =  3 * psfWidth - 2;
    int const paddedHeight = 3 * psfHeight - 2;
    // initialize padded reflected psf image to 0 because only the center is set by computeImage
    afwImage::Image<double> paddedReflPsfImage(paddedWidth, paddedHeight, 0);
    afwImage::BBox reflPsfBBox(afwImage::PointI(psfWidth-1, psfHeight-1), psfWidth, psfHeight);
    afwImage::Image<double> reflPsfImage(paddedReflPsfImage, reflPsfBBox);
    psfKernel.computeImage(reflPsfImage, true);
    reflectImage(reflPsfImage);
    // no need to initialize padded blurred PSF image because all pixels are set
    afwImage::Image<double> paddedBlurredPsfImage(paddedWidth, paddedHeight);
    afwMath::convolve(paddedBlurredPsfImage, paddedReflPsfImage, psfKernel, _normalizePsf);
    afwImage::BBox blurredPsfBBox(afwImage::PointI(psfKernel.getCtrX(), psfKernel.getCtrY()),
        _blurredPsfImage.getWidth(), _blurredPsfImage.getHeight());
    _blurredPsfImage <<= afwImage::Image<double>(paddedBlurredPsfImage, blurredPsfBBox);
};

/**
 * \brief Compute _blurredEposure = scienceExposure convolved with psfKernel
 *
 * \warning If you want scienceExposure convolved with psfKernel(-r)
 * (the standard Kaiser thing to do) then feed in psfKernel(-r)
 *
 * \ingroup coadd::kaiser
 */
void coaddKaiser::CoaddComponent::computeBlurredExposure(
    ExposureF const &scienceExposure,   ///< science exposure
    lsst::afw::math::Kernel const &psfKernel    ///< PSF kernel
) {
    ExposureCC::MaskedImageT blurredMI = _blurredExposure.getMaskedImage();
    ExposureF::MaskedImageT const scienceMI = scienceExposure.getMaskedImage();
//     scienceExposure.writeFits("scienceExposure");
    afwMath::convolve(blurredMI, scienceMI, psfKernel, _normalizePsf);
//     _blurredExposure.writeFits("blurredExposure");
    if (scienceExposure.hasWcs()) {
        afwImage::Wcs::Ptr scienceWcsPtr = scienceExposure.getWcs();
        _blurredExposure.setWcs(*scienceWcsPtr);
    }
};
