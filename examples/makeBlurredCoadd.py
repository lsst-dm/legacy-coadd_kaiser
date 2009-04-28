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

DefSaveImages = True

def makeBlankTemplateExposure(fromExposure):
    """Generate a blank template from a maskedImage
    
    The template will have:
    - Four times as many pixels (twice the resolution in x and y)
    - The same center RA/Dec
    - 2x the on-sky resolution in RA and Dec
    - North up, east to the right
    
    Note that the amount of overlap between the template and maskedImage will vary significantly
    based on how maskedImage is rotated.
    
    Warning:
    - The maskedImage must use FK5 J2000 coordinates for its WCS. This is NOT checked.
    """
    fromMaskedImage = fromExposure.getMaskedImage()
    fromShape = numpy.array(fromMaskedImage.getDimensions(), dtype=int)
    fromWcs = fromExposure.getWcs()

    templateShape = fromShape * 2
    templateMaskedImage = afwImage.MaskedImageD(templateShape[0], templateShape[1])
    templateMaskedImage.set((0,0,0))
    
    # make tangent-plane projection WCS for the template
    fromCtr = (fromShape - 1) / 2
    fromCtrPt = afwImage.PointD(*fromCtr)
    raDecCtr = numpy.array(fromWcs.xyToRaDec(fromCtrPt))
    templateDegPerPix = math.sqrt(fromWcs.pixArea(fromCtrPt))
    
    templateMetadata = dafBase.PropertySet()
    templateMetadata.add("EPOCH", 2000.0)
    templateMetadata.add("EQUINOX", 2000.0)
    templateMetadata.add("RADECSYS", "FK5")
    templateMetadata.add("CTYPE1", "RA---TAN")
    templateMetadata.add("CTYPE2", "DEC--TAN")
    templateMetadata.add("CRPIX1", templateShape[0]/2)
    templateMetadata.add("CRPIX2", templateShape[1]/2)
    templateMetadata.add("CRVAL1", raDecCtr[0])
    templateMetadata.add("CRVAL2", raDecCtr[1])
    templateMetadata.add("CD1_1", templateDegPerPix)
    templateMetadata.add("CD1_2", 0.0)
    templateMetadata.add("CD2_1", 0.0)
    templateMetadata.add("CD2_2", templateDegPerPix)
    templateWcs = afwImage.Wcs(templateMetadata)
    templateExposure = afwImage.ExposureD(templateMaskedImage, templateWcs)
    return templateExposure


def subtractBackground(exposure, doDisplay = False):
    """Subtract the background from an Exposure
    
    Returns the background object returned by afwMath.makeBackground.
    """
    maskedImage = exposure.getMaskedImage()
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    bkgControl.setNxSample(max(2, int(maskedImage.getWidth()/256) + 1))
    bkgControl.setNySample(max(2, int(maskedImage.getHeight()/256) + 1))
    bkgControl.sctrl.setNumSigmaClip(3)
    bkgControl.sctrl.setNumIter(3)
    
    im = maskedImage.getImage()
    bkgObj = afwMath.makeBackground(im, bkgControl)
    im -= bkgObj.getImageF()
    if doDisplay:
        ds9.mtv(maskedImage)
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

where indata is a file containing a list of:
    pathToExposure pathToPsf
where:
    - pathToExposure is the path to an Exposure (without the final _img.fits)
    - pathToPsf is the path to a Psf model persisted as an XML file
Empty lines and lines that start with # are ignored.

The policy controlling the parameters is makeBlurredCoadd_policy.paf
"""
    if len(sys.argv) != 3:
        print helpStr
        sys.exit(0)
    
    outname = sys.argv[1]
    if os.path.exists(outname):
        print "Coadd file %s already exists" % (outname,)
        print helpStr
        sys.exit(1)
    
    indata = sys.argv[2]

    saveImages = DefSaveImages # could make this a command-line option

    makeBlurredCoaddPolicyPath = DefPolicyPath
    makeBlurredCoaddPolicy = pexPolicy.Policy.createPolicy(makeBlurredCoaddPolicyPath)
    detectSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("detectSourcesPolicy")
    psfPolicy = detectSourcesPolicy.getPolicy("psfPolicy")
    measureSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("measureSourcesPolicy")
    fitPsfPolicy = makeBlurredCoaddPolicy.getPolicy("fitPsfPolicy")

    # parse indata
    ImageSuffix = "_img.fits"
    imageDataList = []
    templateExposure = None
    templateMaskedImage = None
    templateWcs = None
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
            maskedImage = exposure.getMaskedImage()
            psfModel = unpersistPsf(psfPath)
            
            if not templateExposure:
                templateExposure = makeBlankTemplateExposure(exposure)
                templateMaskedImage = templateExposure.getMaskedImage()
                templateWcs = templateExposure.getWcs()
            
            print "  Fit and subtract background"
            subtractBackground(exposure)
            if saveImages:
                exposure.writeFits("bgSubtracted%s" % (fileName,))
            
            print "  Compute coadd component"
            coaddComponent = coaddKaiser.CoaddComponent(exposure, psfModel.getKernel())

            print "  Divide exposure by sigma squared = %s" % (coaddComponent.getSigmaSq(),)
            blurredExposure = coaddComponent.getBlurredExposure()
            blurredMaskedImage = blurredExposure.getMaskedImage()
            sigmaSq = coaddComponent.getSigmaSq()
            if saveImages:
                blurredExposure.writeFits("blurred%s" % (fileName,))
            blurredMaskedImage /= sigmaSq
            if saveImages:
                blurredExposure.writeFits("scaledBlurred%s" % (fileName,))
            
            print "  Remap blurred exposure to match template WCS"
            remappedBlurredMaskedImage = afwImage.MaskedImageD(
                templateExposure.getWidth(), templateExposure.getHeight())
            remappedBlurredExposure = afwImage.ExposureD(remappedBlurredMaskedImage, templateWcs)
            if saveImages:
                remappedBlurredExposure.writeFits("remappedBlurred%s" % (fileName,))
            nGoodPix = afwMath.warpExposure(remappedBlurredExposure, blurredExposure,
                afwMath.LanczosWarpingKernel(3))
            nPix = templateExposure.getWidth() * templateExposure.getHeight()
            print "  Remapped image has %d good pixels (%0.0f %%)" % (nGoodPix, 100 * nGoodPix / float(nPix))

            print "  Add remapped blurred exposure to template and save updated template exposure"
            templateMaskedImage += remappedBlurredExposure.getMaskedImage()
            templateMaskedImage.writeFits(outname)
