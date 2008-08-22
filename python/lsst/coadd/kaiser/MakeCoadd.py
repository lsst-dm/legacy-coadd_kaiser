"""
@brief Create Kaiser Coadd

@file
"""
import lsst.coadd.kaiser as coaddKaiser

def makeBlurredCoadd(
    outputExposure,
    scienceExposureList,
    psfFunction,
):
    """Generate a Kaiser coadd
    
    The result is a set of science images convolved
    with their own PSF and then summed
    
    Adds a set of science images convolved with their own PSF
    
    In/Out:
    - outputExposure    coadd
                on input contains desired WCS and size
                on output contains coadd
    
    Inputs:
    - scienceExposureList   sequence of science exposures
    - psfFunction     model for PSF
    
    Returns:
    - ? may want to return science images wcs-mapped
    and convolved with their PSFs, plus associated
    PSF convolved with PSF. Either return this info
    or save it to files to speed up reruns and aid debugging.
    
    Note: in the long run we may want to first wcsRemap each
    science exposure into the WCS of the output image
    (but offset in an Exposure that is large enough and
    centered in the right spot so as to just include
    all the original pixels).
    
    However, that creates very large images (it more than
    quadruples the # of pixels for each science image--
    4x due to scale change and more due to adding a frame of mostly
    garbage pixles to handle rotation). So it is likely to slow things down.
    
    So for now, just try it the "easy" and more compact way. It's likely
    to work except assuming a constant PSF is not smart if the image
    is from a distorted part of the focal plane.
    """
    for scienceExposure in scienceExposureList:
        # compute sigma squared = median of variance
        # (note: might want to compare to value after wcsMatch)
        begPtr = 0 # what?
        endPtr = 0 # what?
        sigSq = coaddKaiser.medianBinapprox(begPtr, endPtr)

        # fit PSF
        # not spatially varying since we don't know how to handle
        # that yet in the Kaiser method
        psfKernel = afwMath.AnalyticKernelD(psfFunction)
        
        # create blurred science subimage =
        # convolve science subimage with PSF kernel

        # WCS-map blurred exposure onto outputExposure
        # we don't really need the remapped pixels to stand alone,
        # but instead would be happy to just get each pixel in turn
        # so we could add it to the output exposure.
        # Still...wcsMatch doesn't work that way at this time.
        
        # add convolved science subimage to coadd
        # weighted by sigma squared
        
        # compute blurred PSF = PSF convolved with itself;
        # 0-pad one of the PSFs first.
        
        # persist or otherwise record blurred science subimage
        # and blurred PSF.