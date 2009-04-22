#!/usr/bin/env python
from __future__ import with_statement
"""
Test:
on my Unix box:
python examples/makeBlurredCoadd.py testTeamplate examples/small.txt /net/scratch1/rowen/cfhtdata/

On my Mac
python examples/makeBlurredCoadd.py testTemplate examples/small.txt /Users/rowen/CFHTData/
"""
import os
import sys
import math
import numpy

from lsst.daf.base import PropertySet
import lsst.pex.logging as pexLog
import lsst.pex.policy as pexPolicy
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDetect
import lsst.afw.display.ds9 as ds9
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.Psf # should be automatically imported, but oh well
import lsst.sdqa as sdqa

RadPerDeg = math.pi / 180.0

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
    fromCtr = (fromShape - 1) / 2
    fromWcs = fromExposure.getWcs()

    templateShape = fromShape * 2
    templateMaskedImage = afwImage.MaskedImageD(templateShape[0], templateShape[1])
    templateMaskedImage.set((0,0,0))
    
    # make tangent-plane projection WCS for the template
    raDecCtr = numpy.array(fromWcs.xyToRaDec(fromCtr[0], fromCtr[1]))
    raDecUp = numpy.array(fromWcs.xyToRaDec(fromCtr[0], fromCtr[1]+1))
    dRaDec = raDecUp - raDecCtr
    dRaDec[0] *= math.cos(raDecCtr[0] * RadPerDeg)
    fromDegPerPix = math.sqrt(dRaDec[0]**2 + dRaDec[1]**2)
    templateDegPerPix = fromDegPerPix / 2.0
    
    templateMetadata = PropertySet()
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


def detectSources(maskedImage, detectSourcesPolicy, doDisplay = False):
    """Measure sources on an Exposure.
    
    This really should be a subroutine in meas_algorithms and I submitted PR #743 requesting that.
    Meanwhile this code is taken from meas_pipelines SourceDetectionStage.
    
    Inputs:
    - maskedImage: masked image on which to measure sources
    - detectSourcesPolicy: a Policy containing elements (with suggested values):
        minPixels:1 
        thresholdValue: 3
        thresholdType: "stdev"
        psfPolicy: {
            algorithm: "DoubleGaussian"
            width = 15
            height = 15
            #5*/(2*sqrt(2*log(2)))
            parameter: 3.22195985
        }
    - doDisplay: True to display diagnostic information on ds9
    
    psfPolicy is used to smooth the image before detecting sources.
    """
    # parse policy
    psfPolicy = detectSourcesPolicy.getPolicy("psfPolicy")
    psf = makePsf(psfPolicy)

    minPixels = detectSourcesPolicy.get("minPixels")
    thresholdValue = detectSourcesPolicy.get("thresholdValue")
    thresholdType = detectSourcesPolicy.get("thresholdType")
    thresholdObj = afwDetect.createThreshold(thresholdValue, thresholdType, True)    
    
    smoothedMaskedImage = maskedImage.Factory(maskedImage.getDimensions())
    smoothedMaskedImage.setXY0(maskedImage.getXY0())

    if doDisplay:
        ds9.mtv(smoothedMaskedImage)
        
    # Smooth the Image
    psf.convolve(smoothedMaskedImage, maskedImage, True, smoothedMaskedImage.getMask().getMaskPlane("EDGE"))

    # Only search psf-smooth part of frame
    llc = afwImage.PointI(psf.getKernel().getWidth()/2,  psf.getKernel().getHeight()/2)
    urc = afwImage.PointI(smoothedMaskedImage.getWidth() - 1, smoothedMaskedImage.getHeight() - 1)
    urc -= llc
    bbox = afwImage.BBox(llc, urc)
    middle = smoothedMaskedImage.Factory(smoothedMaskedImage, bbox)
   
    detectionSet = afwDetect.makeDetectionSet(middle, thresholdObj, "DETECTED", minPixels)
    # detectionSet only searched the middle but it belongs to the entire MaskedImage
    detectionSet.setRegion(afwImage.BBox(afwImage.PointI(maskedImage.getX0(), maskedImage.getY0()),
                                       maskedImage.getWidth(), maskedImage.getHeight()));

    # Grow the detections into the edge by at least one pixel so that it sees the EDGE bit
    grow, isotropic = 1, False
    detectionSet = afwDetect.DetectionSetF(detectionSet, grow, isotropic)
    detectionSet.setMask(maskedImage.getMask(), "DETECTED")

    return detectionSet

def makePsf(psfPolicy):
    params = []        
    params.append(psfPolicy.getString("algorithm"))
    params.append(psfPolicy.getInt("width"))
    params.append(psfPolicy.getInt("height"))
    if psfPolicy.exists("parameter"):
        params += psfPolicy.getDoubleArray("parameter")
    
    return measAlg.createPSF(*params)


def measureSources(exposure, detectionSet, measureSourcesPolicy, doDisplay=False):
    """Measure the characteristics of sources.

    This really should be a subroutine in meas_algorithms and I submitted PR #743 requesting that.
    Meanwhile this code is taken from meas_pipelines SourceMeasurementStage.
    
    measureSourcesPolicy: { 
        centroidAlgorithm: "SDSS"
        shapeAlgorithm: "SDSS"
        photometryAlgorithm: "NAIVE"
        apRadius: 3.0
        psfPolicy: {
            algorithm: "DoubleGaussian"
            width = 15
            height = 15
            #5*/(2*sqrt(2*log(2)))
            parameter: 3.22195985
        }
    }
    """
    psfPolicy = measureSourcesPolicy.getPolicy("psfPolicy")
    psf = makePsf(psfPolicy)
    measureSources = measAlg.makeMeasureSources(exposure, measureSourcesPolicy, psf)
    
    footprints = detectionSet.getFootprints()
    sourceList = afwDetect.SourceSet()
    for i in range(len(footprints)):
        source = afwDetect.Source()
        sourceList.append(source)
    
        source.setId(i)
        source.setFlagForDetection(source.getFlagForDetection() | measAlg.Flags.BINNED1);
    
        try:
            measureSources.apply(source, footprints[i])
        except Exception, e:
            try:
                print e
            except Exception, ee:
                print ee
        
        if source.getFlagForDetection() & measAlg.Flags.EDGE:
            continue
    
        if doDisplay:
            xc, yc = source.getXAstrom() - mi.getX0(), source.getYAstrom() - mi.getY0()
            if False:
                ds9.dot("%.1f %d" % (source.getPsfFlux(), source.getId()), xc, yc+1)
    
            ds9.dot("+", xc, yc, size=1)
            
            if source.getFlagForDetection() & (measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER):
                continue
            if False:               # XPA causes trouble
                Ixx, Ixy, Iyy = source.getIxx(), source.getIxy(), source.getIyy()
                ds9.dot("@:%g,%g,%g" % (Ixx, Ixy, Iyy), xc, yc)
    return sourceList


def fitPsf(exposure, sourceList, fitPsfPolicy):
    """Fit PSF based on a set of measured sources.
    
    fitPsfPolicy must contain:
    - fluxLim
    - sizeCellX
    - sizeCellY
    - nEigenComponents
    - spatialOrder
    - nStarPerCell
    - kernelSize
    - nStarPerCellSpatialFit
    - tolerance
    - reducedChi2ForPsfCandidates
    - nIterForPsf

    Return a fit PSF and a psfCelLSet (whatever that is)
    """
    sdqaRatings = sdqa.SdqaRatingSet()

    psf, psfCellSet = measAlg.Psf.getPsf(
        exposure = exposure,
        sourceList = sourceList,
        moPolicy = fitPsfPolicy,
        sdqaRatings = sdqaRatings,
    )
    return psf, psfCellSet


if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeBlurredCoadd.py coaddfile indata [indir]

where indata is a file containing a list of paths to Exposures (without the final _img.fits);
these paths are relative to indir, if indir is specified.
Empty lines and lines that start with # are ignored.

The policy controlling the parameters is makeBlurredCoadd_policy.paf
"""
    if len(sys.argv) not in (3, 4):
        print helpStr
        sys.exit(0)
    
    outname = sys.argv[1]
    if os.path.exists(outname):
        print "Coadd file %s already exists" % (outname,)
        print helpStr
        sys.exit(1)
    
    indata = sys.argv[2]
    if len(sys.argv) > 3:
        indir = sys.argv[3]
    else:
        indir = "."

    makeBlurredCoaddPolicyPath = "makeBlurredCoadd_policy.paf"
    makeBlurredCoaddPolicy = pexPolicy.Policy.createPolicy(makeBlurredCoaddPolicyPath)
    detectSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("detectSourcesPolicy")
    psfPolicy = detectSourcesPolicy.getPolicy("psfPolicy")
    measureSourcesPolicy = makeBlurredCoaddPolicy.getPolicy("measureSourcesPolicy")
    fitPsfPolicy = makeBlurredCoaddPolicy.getPolicy("fitPsfPolicy")

    # parse indata
    ImageSuffix = "_img.fits"
    imageDataList = []
    templateExposure = None
    with file(indata, "rU") as infile:
        for lineNum, line in enumerate(infile):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            filename = line
            filepath = os.path.join(indir, filename)
            if not os.path.isfile(filepath + ImageSuffix):
                print "Skipping exposure %s; file %s not found" % (filepath, filepath + ImageSuffix)
                continue
            
            print "Processing exposure %s" % (filepath,)
            exposure = afwImage.ExposureF(filepath)
            maskedImage = exposure.getMaskedImage()
            
            if not templateExposure:
                templateExposure = makeBlankTemplateExposure(exposure)
            
            subtractBackground(exposure)
            
            # fit a spatially varying PSF
            detectionSet = detectSources(maskedImage, detectSourcesPolicy, doDisplay = False)
            sourceList = measureSources(exposure, detectionSet, measureSourcesPolicy, doDisplay=False)
            psf, psfCellSet = fitPsf(exposure, sourceList, fitPsfPolicy)

            # compute coadd component
            
            # add coadd component to template
            