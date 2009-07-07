#!/usr/bin/env python
from __future__ import with_statement
"""
This example requires:
- A set of science exposures
- A set PSF models for those exposures, persisted as XML files
- A file containing the paths to each, as:
  image1 psf1
  image2 psf2
  ...
  
Test:
On my Mac
python examples/makeBlurredCoadd.py testTemplate examples/imagesToCoadd.txt
"""
import os
import sys
import math
import numpy

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.Psf # should be automatically imported, but oh well
import lsst.sdqa as sdqa
import lsst.coadd.kaiser as coaddKaiser

BaseDir = os.path.dirname(__file__)
DefPolicyPath = os.path.join(BaseDir, "makeBlurredCoadd_policy.paf")

RadPerDeg = math.pi / 180.0

DefSaveImages = False
# names of mask planes permitted in the coadd
AcceptableMaskPlaneList = ("SAT", "INTRP")

def makeBlankTemplateExposure(fromExposure):
    """Generate a blank coadd from a maskedImage
    
    The coadd will have:
    - Four times as many pixels (twice the resolution in x and y)
    - The same center RA/Dec
    - 2x the on-sky resolution in RA and Dec
    - North up, east to the right
    
    Note that the amount of overlap between the coadd and maskedImage will vary significantly
    based on how maskedImage is rotated.
    
    Warning:
    - The maskedImage must use FK5 J2000 coordinates for its WCS. This is NOT checked.
    """
    fromMaskedImage = fromExposure.getMaskedImage()
    fromShape = numpy.array(fromMaskedImage.getDimensions(), dtype=int)
    fromWcs = fromExposure.getWcs()

    coaddShape = fromShape * 2
    coaddMaskedImage = afwImage.MaskedImageD(coaddShape[0], coaddShape[1])
    coaddMaskedImage.set((0,0,0))
    
    # make tangent-plane projection WCS for the coadd
    fromCtr = (fromShape - 1) / 2
    fromCtrPt = afwImage.PointD(*fromCtr)
    raDecCtr = numpy.array(fromWcs.xyToRaDec(fromCtrPt))
    coaddDegPerPix = math.sqrt(fromWcs.pixArea(fromCtrPt)) / 2.0
    
    coaddMetadata = dafBase.PropertySet()
    coaddMetadata.add("EPOCH", 2000.0)
    coaddMetadata.add("EQUINOX", 2000.0)
    coaddMetadata.add("RADECSYS", "FK5")
    coaddMetadata.add("CTYPE1", "RA---TAN")
    coaddMetadata.add("CTYPE2", "DEC--TAN")
    coaddMetadata.add("CRPIX1", coaddShape[0]/2)
    coaddMetadata.add("CRPIX2", coaddShape[1]/2)
    coaddMetadata.add("CRVAL1", raDecCtr[0])
    coaddMetadata.add("CRVAL2", raDecCtr[1])
    coaddMetadata.add("CD1_1", coaddDegPerPix)
    coaddMetadata.add("CD1_2", 0.0)
    coaddMetadata.add("CD2_1", 0.0)
    coaddMetadata.add("CD2_2", coaddDegPerPix)
    coaddWcs = afwImage.Wcs(coaddMetadata)
    coaddExposure = afwImage.ExposureD(coaddMaskedImage, coaddWcs)
    return coaddExposure


def subtractBackground(maskedImage, doDisplay = False):
    """Subtract the background from a MaskedImage
    
    Note: at present the mask and variance are ignored, but they might used be someday.
    
    Returns the background object returned by afwMath.makeBackground.
    """
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    bkgControl.setNxSample(max(2, int(maskedImage.getWidth()/256) + 1))
    bkgControl.setNySample(max(2, int(maskedImage.getHeight()/256) + 1))
    bkgControl.sctrl.setNumSigmaClip(3)
    bkgControl.sctrl.setNumIter(3)

    image = maskedImage.getImage()
    bkgObj = afwMath.makeBackground(image, bkgControl)
    image -= bkgObj.getImageF()
    if doDisplay:
        ds9.mtv(image)
    return bkgObj

def unpersistPsf(xmlPath):
    """Read a PSF from an XML file"""
    # Set up persistence object
    pol = pexPolicy.Policy()
    persistence = dafPersist.Persistence.getPersistence(pol)
    
    # Where is the file on disk? Make a storage object
    loc = dafPersist.LogicalLocation(xmlPath)
    storageList = dafPersist.StorageList()
    storage= persistence.getRetrieveStorage('XmlStorage', loc)
    storageList.append(storage)
    
    # Capture any associated metadata
    metadata = dafBase.PropertySet()
    
    # Unpersist the object; you need to say which object in storage to grab
    persistable = persistence.retrieve('pcaPsf', storageList, metadata)
    
    # Cast to a PSF model
    psf = measAlg.PSF_swigConvert(persistable)
    return psf
    

if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeBlurredCoadd.py coaddfile indata

where:
- coaddfile is the desired name or path of the output coadd
- indata is a file containing a list of:
  pathToExposure pathToPsf
  where:
  - pathToExposure is the path to an Exposure (without the final _img.fits)
  - pathToPsf is the path to a Psf model persisted as an XML file
  - empty lines and lines that start with # are ignored.

The policy controlling the parameters is makeBlurredCoadd_policy.paf
"""
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    outName = sys.argv[1]
    if os.path.exists(outName + "_img.fits"):
        print "Coadd file %s already exists" % (outName,)
        print helpStr
        sys.exit(1)
    depthOutName = outName + "_depth.fits"
    
    indata = sys.argv[2]
    
    acceptableMask = 0
    for maskPlaneName in AcceptableMaskPlaneList:
        acceptableMask |= 1 << afwImage.MaskU.getMaskPlane(maskPlaneName)
    CoaddMask = 0xFFFF - acceptableMask
    

    saveImages = DefSaveImages # could make this a command-line option

    makeBlurredCoaddPolicyPath = DefPolicyPath
    makeBlurredCoaddPolicy = pexPolicy.Policy.createPolicy(makeBlurredCoaddPolicyPath)
    normalizePsf = makeBlurredCoaddPolicy.get("normalizePsf")
    detectSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("detectSourcesPolicy")
    psfPolicy = detectSourcesPolicy.getPolicy("psfPolicy")
    measureSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("measureSourcesPolicy")
    fitPsfPolicy = makeBlurredCoaddPolicy.getPolicy("fitPsfPolicy")

    # parse indata
    ImageSuffix = "_img.fits"
    imageDataList = []
    coaddExposure = None
    coaddMaskedImage = None
    coaddWcs = None
    depthMap = None
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            filePath, psfPath = line.split()[0:2]
            fileName = os.path.basename(filePath)
            if not os.path.isfile(filePath + ImageSuffix):
                print "Skipping exposure %s; image file %s not found" % (fileName, filePath + ImageSuffix,)
                continue
            if not os.path.isfile(psfPath):
                print "Skipping exposure %s; psf file %s not found" % (fileName, psfPath)
                continue
            
            print "Processing exposure %s" % (filePath,)
            exposure = afwImage.ExposureF(filePath)
#           Uncomment to make a quick small template
#             subStart = afwImage.PointI(0, 0)
#             subExposureBBox = afwImage.BBox(subStart, 300, 300)
#             exposure = afwImage.ExposureF(exposure, subExposureBBox, True)
            maskedImage = exposure.getMaskedImage()
            psfModel = unpersistPsf(psfPath)
            psfKernel = psfModel.getKernel()
            if normalizePsf and afwMath.LinearCombinationKernel.swigConvert(psfKernel) != None:
                # psf kernel is a LinearCombinationKernel; convert to a FixedKernel
                if psfKernel.isSpatiallyVarying():
                    print """Warning: ignoring spatial variation of psf: normalizePsf is True,
but psf kernel is a spatially varying LinearCombinationKernel,
which cannot be normalized until ticket 833 is implemented."""
                psfImage = afwImage.ImageD(psfKernel.getWidth(), psfKernel.getHeight())
                xCtrInd = (psfKernel.getWidth() - 1) / 2.0
                yCtrInd = (psfKernel.getHeight() - 1) / 2.0
#               note: use indexToPosition once a floating point version is available -- ticket #845
#                 xCtrPos = afwImage.indexToPosition(xCtrInd)
#                 yCtrPos = afwImage.indexToPosition(yCtrInd)
                xCtrPos = afwImage.PixelZeroPos + xCtrInd
                yCtrPos = afwImage.PixelZeroPos + yCtrInd
                psfKernel.computeImage(psfImage, True, xCtrPos, yCtrPos)
                psfKernel = afwMath.FixedKernel(psfImage)
            
            if not coaddExposure:
                coaddExposure = makeBlankTemplateExposure(exposure)
                coaddMaskedImage = coaddExposure.getMaskedImage()
                coaddWcs = coaddExposure.getWcs()
                depthMap = afwImage.ImageU(coaddMaskedImage.getDimensions(), 0)
            
            print "  Fit and subtract background"
            subtractBackground(maskedImage)
            if saveImages:
                exposure.writeFits("bgSubtracted%s" % (fileName,))
            
            print "  Compute coadd component"
            coaddComponent = coaddKaiser.CoaddComponent(exposure, psfKernel, normalizePsf)

            print "  Divide exposure by sigma squared = %s" % (coaddComponent.getSigmaSq(),)
            blurredExposure = coaddComponent.getBlurredExposure()
            blurredMaskedImage = blurredExposure.getMaskedImage()
            sigmaSq = coaddComponent.getSigmaSq()
            if saveImages:
                blurredExposure.writeFits("blurred%s" % (fileName,))
            blurredMaskedImage /= sigmaSq
            if saveImages:
                blurredExposure.writeFits("scaledBlurred%s" % (fileName,))
            
            print "  Remap blurred exposure to match coadd WCS"
            remappedBlurredMaskedImage = afwImage.MaskedImageD(
                coaddExposure.getWidth(), coaddExposure.getHeight())
            remappedBlurredExposure = afwImage.ExposureD(remappedBlurredMaskedImage, coaddWcs)
            if saveImages:
                remappedBlurredExposure.writeFits("remappedBlurred%s" % (fileName,))
            nGoodPix = afwMath.warpExposure(remappedBlurredExposure, blurredExposure,
                afwMath.LanczosWarpingKernel(3))
            nPix = coaddExposure.getWidth() * coaddExposure.getHeight()
            print "  Remapped image has %d good pixels (%0.0f %%)" % (nGoodPix, 100 * nGoodPix / float(nPix))

            print "  Add remapped blurred exposure to coadd and save updated coadd exposure"
            coaddKaiser.addToCoadd(coaddMaskedImage, depthMap, remappedBlurredExposure.getMaskedImage(), CoaddMask)
            coaddExposure.writeFits(outName)
            depthMap.writeFits(depthOutName)
    coaddKaiser.setCoaddEdgeBits(coaddMaskedImage.getMask(), depthMap)
    coaddExposure.writeFits(outName)    