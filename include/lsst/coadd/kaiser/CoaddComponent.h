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
    * @brief One component (processed Exposure) of a Kaiser coadd
    */
    class CoaddComponent : public lsst::daf::data::LsstBase {
    public:
        typedef float pixelType;
        typedef lsst::afw::image::Exposure<pixelType, lsst::afw::image::maskPixelType> Exposure;
        typedef lsst::afw::image::MaskedImage<pixelType, lsst::afw::image::maskPixelType> MaskedImage;

        CoaddComponent();
        CoaddComponent(
            Exposure const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
        virtual ~CoaddComponent() {};

        void addToCoadd(Exposure &coadd);
        
        double getSigmaSq() { return _sigmaSq; }

        Exposure getBlurredExposure() { return _blurredExposure; }
        
        lsst::afw::image::Image<pixelType> getBlurredPsfImage() { return _blurredPsfImage; };
        
    private:
        double _sigmaSq;
        Exposure _blurredExposure;
        lsst::afw::image::Image<pixelType> _blurredPsfImage;
        
        void computeSigmaSq(
            Exposure const &scienceExposure
        );
        
        void computeBlurredPsf(
            lsst::afw::math::Kernel const &psfKernel
        );

        void computeBlurredExposure(
            Exposure const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
    };

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H)
