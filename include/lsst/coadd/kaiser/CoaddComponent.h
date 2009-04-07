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
     * @brief One component (processed ExposureD) of a Kaiser coadd
     *
     * @ingroup coadd::kaiser
     */
    class CoaddComponent : public lsst::daf::data::LsstBase {
    public:
        typedef double pixelType; // pixel type for blurred science exposure
        typedef lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel> ExposureF;
        typedef lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel> MaskedImageF;
        typedef lsst::afw::image::Exposure<pixelType, lsst::afw::image::MaskPixel> ExposureD;
        typedef lsst::afw::image::MaskedImage<pixelType, lsst::afw::image::MaskPixel> MaskedImageD;

        CoaddComponent(
            ExposureF const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
        virtual ~CoaddComponent() {};

        void addToCoadd(ExposureD &coadd);
        
        double getSigmaSq() { return _sigmaSq; }

        ExposureD getBlurredExposure() { return _blurredExposure; }
        
        lsst::afw::image::Image<double> getBlurredPsfImage() { return _blurredPsfImage; };
        
    private:
        double _sigmaSq;
        ExposureD _blurredExposure;
        lsst::afw::image::Image<double> _blurredPsfImage;
        
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

    void reflectImage(lsst::afw::image::Image<CoaddComponent::pixelType> &image);

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H)
