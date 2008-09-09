"""
@brief Create Kaiser Coadd

@file
"""
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import lsst.coadd.kaiser as coaddKaiser

def makeBlurredCoadd(
    outputExposure,
    scienceExposureList,
    psfKernelList,
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
    - psfKernelList      list of psf kernels
    
    At the moment it discards the CoadComponents that it generates.
    We'll probably need to save these somehow, but I'm afraid if we return them
    we'll run out of memory.
    """
    if len(scienceExposureList) != len(psfKernelList):
        raise RuntimeError("len(scienceExposureList) != len(psfKernelList)")

    RemapKernelSize = 5
    outputMaskedImage = outputExposure.getMaskedImage()
    tempMaskedImage = afwImage.MaskedImageD(outputMaskedImage.getCols(), outputMaskedImage.getRows())
    tempExposure = afwImage.ExposureD(tempMaskedImage, outputExposure.getWCS())

    for scienceExposure, psfKernel in zip(scienceExposureList, psfKernelList):
        coaddComp = coaddKaiser.CoaddComponent(scienceExposure, psfKernel)

        # WCS-map blurred exposure onto outputExposure
        ipDiffim.wcsMatch(tempExposure, coaddComp.getBlurredExposure(), "lanczos", RemapKernelSize, RemapKernelSize)
        tempMaskedImage = tempExposure.getMaskedImage()
        weightingFactor = 1.0 / coaddComp.getSigmaSq()
        tempMaskedImage *= weightingFactor
        
        # add convolved science subimage to coadd
        # weighted by sigma squared
        # do NOT add masked pixels from blurredExposure
        badPixelMask = 0xFFFF
        coaddKaiser.addToMaskedImage(outputMaskedImage, tempMaskedImage, badPixelMask)
