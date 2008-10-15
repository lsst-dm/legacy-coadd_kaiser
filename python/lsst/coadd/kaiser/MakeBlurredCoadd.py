"""
@brief Create Kaiser Coadd

@file
"""
__all__ = ["makeBlurredCoadd"]

import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.ip.diffim as ipDiffim
import kaiserLib

def makeBlurredCoadd(
    outputExposure,
    imageDataList,
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
    - imageDataList     a list of objects that contain at least these attributes:
        - filepath      full path to image file, excluding "_img.fits"
        - sky           sky level to subtract from the image
        - psfKernel     a PSF kernel model
    
    At the moment it discards the CoadComponents that it generates.
    We'll probably need to save these somehow, but I'm afraid if we return them
    we'll run out of memory.
    
    @todo use a more general mechanism to unpersist the science Exposures
    (but try to avoid inputting all science Exposures at once to save memory).
    """
    RemapKernelSize = 5
    outputMaskedImage = outputExposure.getMaskedImage()
    tempMaskedImage = afwImage.MaskedImageD(outputMaskedImage.getCols(), outputMaskedImage.getRows())
    tempExposure = afwImage.ExposureD(tempMaskedImage, outputExposure.getWcs())

    for imageData in imageDataList:
        pexLog.Trace('lsst.coadd.kaiser.makeBlurredCoadd', 4, "process image %s" % (imageData.filepath,))
        scienceExposure = afwImage.ExposureF()
        scienceExposure.readFits(imageData.filepath)
        if imageData.sky != 0:
            scienceMaskedImage = scienceExposure.getMaskedImage()
            scienceImage = scienceMaskedImage.getImage()
            scienceImage -= imageData.sky
    
        pexLog.Trace('lsst.coadd.kaiser.makeBlurredCoadd', 4, "compute coadd component")
        coaddComp = kaiserLib.CoaddComponent(scienceExposure, imageData.psfKernel)

        # WCS-map blurred exposure onto outputExposure
        pexLog.Trace('lsst.coadd.kaiser.makeBlurredCoadd', 4, "wcs match coadd component")
        ipDiffim.wcsMatch(tempExposure, coaddComp.getBlurredExposure(), "lanczos", RemapKernelSize, RemapKernelSize)
        
        # add convolved science subimage to coadd
        # weighted by sigma squared
        # do NOT add masked pixels from blurredExposure
        tempMaskedImage = tempExposure.getMaskedImage()
        weightingFactor = 1.0 / coaddComp.getSigmaSq()
        tempMaskedImage *= weightingFactor
        badPixelMask = 0xFFFF
        pexLog.Trace('lsst.coadd.kaiser.makeBlurredCoadd', 4, "add coadd component to template")
        kaiserLib.addToMaskedImage(outputMaskedImage, tempMaskedImage, badPixelMask)
