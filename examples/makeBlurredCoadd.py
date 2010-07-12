#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

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
import lsst.coadd.utils as coaddUtils
import lsst.coadd.kaiser as coaddKaiser

BaseDir = os.path.dirname(__file__)
DefPolicyPath = os.path.join(BaseDir, "makeBlurredCoadd_policy.paf")

RadPerDeg = math.pi / 180.0

DefSaveImages = False
# names of mask planes permitted in the coadd
AcceptableMaskPlaneList = ("SAT", "INTRP")

BackgroundCells = 256

def subtractBackground(maskedImage, doDisplay = False):
    """Subtract the background from a MaskedImage
    
    Note: at present the mask and variance are ignored, but they might used be someday.
    
    Returns the background object returned by afwMath.makeBackground.
    """
    if doDisplay:
        ds9.mtv(maskedImage)
    bkgControl = afwMath.BackgroundControl(afwMath.NATURAL_SPLINE)
    bkgControl.setNxSample(int(maskedImage.getWidth() // BackgroundCells) + 1)
    bkgControl.setNySample(int(maskedImage.getHeight() // BackgroundCells) + 1)
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
    

    saveImages = DefSaveImages # could make this a command-line option

    makeBlurredCoaddPolicyPath = DefPolicyPath
    makeBlurredCoaddPolicy = pexPolicy.Policy.createPolicy(makeBlurredCoaddPolicyPath)
    normalizePsf = makeBlurredCoaddPolicy.get("normalizePsf")
    resolutionFactor = policy.get("resolutionFactor")
    allowedMaskPlanes = policy.get("allowedMaskPlanes")
    detectSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("detectSourcesPolicy")
    psfPolicy = detectSourcesPolicy.getPolicy("psfPolicy")
    measureSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("measureSourcesPolicy")
    fitPsfPolicy = makeBlurredCoaddPolicy.getPolicy("fitPsfPolicy")

    acceptableMask = 0
    for maskPlaneName in allowedMaskPlanes:
        acceptableMask |= 1 << afwImage.MaskU.getMaskPlane(maskPlaneName)
    coaddMask = 0xFFFF - acceptableMask
    

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
                coaddExposure = coaddUtils.makeBlankCoadd(exposure, resolutionFactor=resolutionFactor)
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
            coaddUtils.addToCoadd(coaddMaskedImage, depthMap, remappedBlurredExposure.getMaskedImage(),
                coaddMask)
            coaddExposure.writeFits(outName)
            depthMap.writeFits(depthOutName)
    coaddUtils.setCoaddEdgeBits(coaddMaskedImage.getMask(), depthMap)
    coaddExposure.writeFits(outName)    