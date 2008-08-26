// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_PSFKERNEL_H
#define LSST_COADD_KAISER_PSFKERNEL_H
/**
* @brief define double gaussian PSF kernel
*
* @file
*
* @author Russell Owen
*
* @todo move to lsst/afw/math/FunctionLibrary.h or somewhere more relevant
* (e.g. detection, but that has to be coordinated since it already has a
* double gaussian that can't be used to make a kernel).
*/
#include "lsst/afw/math/Function.h"

namespace lsst {
namespace coadd {
namespace kaiser {
    
    /**
     * @brief double Guassian (sum of two Gaussians)
     *
     * Intended for use as a PSF model: the main Gaussian represents the core
     * and the second Gaussian represents the wings.
     *
     * f(x,y) = e^(-r^2 / sigma1^2) + b * e^(-r^2 / sigma2^2)
     * where r^2 = x^2 + y^2
     */
    template<typename ReturnT>
    class DoubleGaussianFunction2: public lsst::afw::math::Function2<ReturnT> {
    public:
        typedef lsst::afw::math::Function2<ReturnT> Function2;
        typedef typename Function2::Ptr Function2Ptr;

        /**
         * @brief Construct a Gaussian function with specified x and y sigma
         */
        explicit DoubleGaussianFunction2(
            double sigma1,      ///< sigma of main Gaussian (which has amplitude 1)
            double sigma2 = 0,  ///< sigma of second Gaussian
            double ampl2 = 0)   ///< amplitude of second Gaussian
        : 
            Function2(2)
        {
            this->_params[0] = sigma1;
            this->_params[1] = sigma2;
            this->_params[2] = ampl2;
        }
        
        virtual ~DoubleGaussianFunction2() {};
        
        virtual Function2Ptr copy() const {
            return Function2Ptr(new DoubleGaussianFunction2(this->_params[0], this->_params[1]));
        }
        
        virtual ReturnT operator() (double x, double y) const {
            double radSq = (x * x) + (y * y);
            return std::exp(-radSq / (2.0 * this->_params[0] * this->_params[0]))
                + this->_params[2] * std::exp(-radSq / (2.0 * this->_params[1] * this->_params[1]));
        }
        
        virtual std::string toString(void) const {
            std::ostringstream os;
            os << "DoubleGaussianFunction2 []: ";
            os << Function2::toString();
            return os.str();
        };
    };

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_PSFKERNEL_H)
