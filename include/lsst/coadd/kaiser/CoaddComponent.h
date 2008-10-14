// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H
#define LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H
/**
* @brief Component of Kaiser coadd
*
* @file
*
* @author Russell Owen
*/
#include <vector>

#include "lsst/daf/data/LsstBase.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

namespace lsst {
namespace coadd {
namespace kaiser {
    
    /**
     * @brief One component (processed lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType>) of a Kaiser coadd
     *
     * @ingroup coadd::kaiser
     *
     * @todo use typedefs for arguments instead of writing out everything
     * but this requires swig 1.3.36+1 and I'm stuck at 1.3.34 for now
     */
    class CoaddComponent : public lsst::daf::data::LsstBase {
    public:
        typedef double pixelType; // pixel type for blurred science exposure
        typedef lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> ExposureF;
        typedef lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> MaskedImageF;
        typedef lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType> ExposureD;
        typedef lsst::afw::image::MaskedImage<pixelType, lsst::afw::image::maskPixelType> MaskedImageD;

        CoaddComponent(
            lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
        virtual ~CoaddComponent() {};

        void addToCoadd(lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType> &coadd);
        
        double getSigmaSq() { return _sigmaSq; }

        lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType> getBlurredExposure() { return _blurredExposure; }
        
        lsst::afw::image::Image<double> getBlurredPsfImage() { return _blurredPsfImage; };
        
    private:
        double _sigmaSq;
        lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType> _blurredExposure;
        lsst::afw::image::Image<double> _blurredPsfImage;
        
        void computeSigmaSq(
            lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &scienceExposure
        );
        
        void computeBlurredPsf(
            lsst::afw::math::Kernel const &psfKernel
        );

        void computeBlurredExposure(
            lsst::afw::image::Exposure<float, lsst::afw::image::maskPixelType> const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
    };

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H)
