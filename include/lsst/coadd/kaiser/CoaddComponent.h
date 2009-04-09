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
#include "lsst/daf/data/LsstBase.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

namespace lsst {
namespace coadd {
namespace kaiser {
    
    /**
     * @brief One component (processed Exposure) of a Kaiser coadd
     *
     * @ingroup coadd::kaiser
     */
    class CoaddComponent : public lsst::daf::data::LsstBase {
    public:
        typedef double pixelType; // pixel type for blurred science exposure
        typedef lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> ExposureF;
        typedef lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> MaskedImageF;
        typedef lsst::afw::image::Exposure<pixelType, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> ExposureCC;
        typedef lsst::afw::image::Image<pixelType> ImageCC;

        CoaddComponent(
            ExposureF const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
        virtual ~CoaddComponent() {};

        void addToCoadd(ExposureCC const &coadd);
        
        double getSigmaSq() const { return _sigmaSq; }

        ExposureCC getBlurredExposure() const { return _blurredExposure; }
        
        ImageCC getBlurredPsfImage() const { return _blurredPsfImage; };
        
    private:
        double _sigmaSq;
        ExposureCC _blurredExposure;
        ImageCC _blurredPsfImage;
        
        void computeSigmaSq(
            ExposureF const &scienceExposure
        );
        
        void computeBlurredPsf(
            lsst::afw::math::Kernel const &psfKernel
        );

        void computeBlurredExposure(
            ExposureF const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
    };

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H)
