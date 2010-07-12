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
 
#ifndef LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H
#define LSST_COADD_KAISER_MAKEBLURREDTEMPLATE_H
/**
* @brief Component of Kaiser coadd
*
* @file
*
* @author Russell Owen
*/
#include "lsst/daf/base/Citizen.h"
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
    class CoaddComponent : public lsst::daf::base::Citizen {
    public:
        typedef double pixelType; // pixel type for blurred science exposure
        typedef lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> ExposureF;
        typedef lsst::afw::image::MaskedImage<float, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> MaskedImageF;
        typedef lsst::afw::image::Exposure<pixelType, lsst::afw::image::MaskPixel,
            lsst::afw::image::VariancePixel> ExposureCC;
        typedef lsst::afw::image::Image<pixelType> ImageCC;

        explicit CoaddComponent(
            ExposureF const &scienceExposure,
            lsst::afw::math::Kernel const &psfKernel,
            bool normalizePsf = true
        );
        virtual ~CoaddComponent() {};

        double getSigmaSq() const { return _sigmaSq; }

        ExposureCC getBlurredExposure() const { return _blurredExposure; }
        
        ImageCC getBlurredPsfImage() const { return _blurredPsfImage; };
        
    private:
        double _sigmaSq;
        ExposureCC _blurredExposure;
        ImageCC _blurredPsfImage;
        bool _normalizePsf;
        
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
