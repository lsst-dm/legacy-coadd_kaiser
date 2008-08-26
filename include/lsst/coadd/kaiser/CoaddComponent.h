// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_makeBlurredTemplate_H
#define LSST_COADD_KAISER_makeBlurredTemplate_H
/**
* @brief Component of Kaiser coadd
*
* @file
*
* @author Russell Owen
*/
#include <vector>

#include "lsst/afw/image.h"
#include "lsst/afw/math.h"

namespace lsst {
namespace coadd {
namespace kaiser {
    
    /**
    * @brief One component (processed Exposure) of a Kaiser coadd
    */
    class CoaddComponent {
    public:
        typedef float pixelType;
        typedef lsst::afw::image::MaskedImage<pixelType, lsst::afw::image::maskPixelType> Exposure;

        CoaddComponent();
        CoaddComponent(
            Exposure const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel
        );
        
        void addToCoadd(Exposure &coadd);
        
        double getSigmaSq() { return _sigmaSq; }

        Exposure getBlurredExposure() { return _blurredExposure; }
        
        lsst::afw::image::Image<pixelType> getBlurredPsf() { return _blurredPsf; };
        
    private:
        double _sigmaSq;
        Exposure _blurredExposure;
        lsst::afw::image::Image<pixelType> _blurredPsf;
        
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

#endif // !defined(LSST_COADD_KAISER_makeBlurredTemplate_H)
