#!/usr/bin/env python
"""
Test:
on my Unix box:
python examples/makeBlurredCoadd.py testTeamplate examples/small.txt /net/scratch1/rowen/cfhtdata/

On my Mac
python examples/makeBlurredCoadd.py testTemplate examples/small.txt /Users/rowen/CFHTData/

Both fail with:
Ready to process the following files:
./729994p_12; med=464.6; fwhm=4.0; kernelSize=16
./729995p_12; med=469.2; fwhm=4.0; kernelSize=16
./729996p_12; med=474.9; fwhm=4.0; kernelSize=16
./729997p_12; med=486.4; fwhm=4.2; kernelSize=17
Generating the blank template from the first file
Building template
    process image ./729994p_12
    compute coadd component
Traceback (most recent call last):
  File "/Users/rowen/LSST/code/coadd_kaiser-trunk/examples/makeBlurredCoadd.py", line 164, in <module>
    coaddKaiser.makeBlurredCoadd(templateExposure, imageDataList)
  File "/Users/rowen/LSST/code/coadd_kaiser-trunk/python/lsst/coadd/kaiser/MakeBlurredCoadd.py", line 55, in makeBlurredCoadd
    coaddComp = kaiserLib.CoaddComponent(scienceExposure, imageData.psfKernel)
  File "/Users/rowen/LSST/code/coadd_kaiser-trunk/python/lsst/coadd/kaiser/kaiserLib.py", line 660, in __init__
    this = _kaiserLib.new_CoaddComponent(*args)
TypeError: in method 'new_CoaddComponent', argument 1 of type 'lsst::coadd::kaiser::CoaddComponent::ExposureF const &'

This appears to NOT be the SWIG bug fixed in 1.3.36+1. Damn.
"""
import os
import sys
import math
import numpy

from lsst.daf.base import DataProperty
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.coadd.kaiser as coaddKaiser

RadPerDeg = math.pi / 180.0

class ImageData(object):
    """Data about a science image
    
    Inputs:
    - filepath      full path to masked image, without the "_img.fits"
    - sky           sky background level (assumed constant over the whole image)
    - fwhm          FWHM of stars on image (assumed constant over the whole image)
    - kernelSize    size of PSF kernel (#rows = #cols); if None then set to 4 x fwhm

    Computed attributes:
    - psfFunction   PSF model as a double Gaussian, the sum of:
        a Gaussian of width FWHM for the central peak
        a Gaussian of width 2*FWHM and 1/10 amplitude for the wings
    - psfKernel     PSF kernel of size kernelSize using psfFunction
    """
    def __init__(self, filepath, sky, fwhm, kernelSize=None):
        self.filepath = filepath
        self.sky = float(sky)
        self.fwhm = float(fwhm)
        sigma = self.fwhm / 2.35
        if kernelSize == None:
            kernelSize = self.fwhm * 4.0
        self.kernelSize = int(round(kernelSize))
        self.psfFunction = afwMath.DoubleGaussianFunction2D(sigma, sigma * 2.0, 0.1)
        self.psfKernel = afwMath.AnalyticKernel(self.psfFunction, self.kernelSize, self.kernelSize)
    
    def __str__(self):
        return "%s; med=%0.1f; fwhm=%0.1f; kernelSize=%s" % (self.filepath, self.sky, self.fwhm, self.kernelSize)

def arrayFromCoord2(coord2):
    return numpy.array([coord2.x(), coord2.y()])

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
    fromShape = numpy.array([fromMaskedImage.getCols(), fromMaskedImage.getRows()], dtype=int)
    fromCtr = (fromShape - 1) / 2
    fromWcs = fromExposure.getWcs()

    templateShape = fromShape * 2
    templateMaskedImage = afwImage.MaskedImageF(templateShape[0], templateShape[1])
    
    # make tangent-plane projection WCS for the template
    raDecCtr = arrayFromCoord2(fromWcs.colRowToRaDec(fromCtr[0], fromCtr[1]))
    raDecUp = arrayFromCoord2(fromWcs.colRowToRaDec(fromCtr[0], fromCtr[1]+1))
    dRaDec = raDecUp - raDecCtr
    dRaDec[0] *= math.cos(raDecCtr[0] * RadPerDeg)
    fromDegPerPix = math.sqrt(dRaDec[0]**2 + dRaDec[1]**2)
    templateDegPerPix = fromDegPerPix / 2.0
    
    templateMetaData = DataProperty.createPropertyNode("FitsMetaData")
    templateMetaData.addProperty(DataProperty("EPOCH", 2000.0))
    templateMetaData.addProperty(DataProperty("EQUINOX", 2000.0))
    templateMetaData.addProperty(DataProperty("RADECSYS", "FK5"))
    templateMetaData.addProperty(DataProperty("CTYPE1", "RA---TAN"))
    templateMetaData.addProperty(DataProperty("CTYPE2", "DEC--TAN"))
    templateMetaData.addProperty(DataProperty("CRPIX1", templateShape[0]/2))
    templateMetaData.addProperty(DataProperty("CRPIX2", templateShape[1]/2))
    templateMetaData.addProperty(DataProperty("CRVAL1", raDecCtr[0]))
    templateMetaData.addProperty(DataProperty("CRVAL2", raDecCtr[1]))
    templateMetaData.addProperty(DataProperty("CD1_1", templateDegPerPix))
    templateMetaData.addProperty(DataProperty("CD1_2", 0.0))
    templateMetaData.addProperty(DataProperty("CD2_1", 0.0))
    templateMetaData.addProperty(DataProperty("CD2_2", templateDegPerPix))
    templateWcs = afwImage.Wcs(templateMetaData)
    templateExposure = afwImage.ExposureF(templateMaskedImage, templateWcs)
    return templateExposure

    
if __name__ == "__main__":
    pexLog.Trace.setVerbosity('lsst.coadd', 5)
    helpStr = """Usage: makeBlurredCoadd.py coaddfile indata [indir]

where indata is a file containing lines of this form:
imagename  median  stdDev  #stars  FWHM
(though only imagename, median, stdDev and FWHM are used)
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

    # parse indata
    imageDataList = []
    infile = file(indata, "rU")
    for line in infile:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        lineData = line.split()
        if len(lineData) < 5:
            print "Could not parse line: %r; skipping it" % (line,)
            continue
        filename = lineData[0]
        if not filename.endswith("_img.fits"):
            print "Skipping %s; not an image file" % (filename,)
            continue
        filepath = os.path.join(indir, filename)
        if not os.path.isfile(filepath):
            print "Skipping image %s; file not found" % (filepath,)
            continue
        filepath = filepath[0:-9]
        sky = lineData[1]
        fwhm = lineData[4]
        imageDataList.append(ImageData(filepath, sky, fwhm))

    print "Ready to process the following files:"
    for imageData in imageDataList:
        print imageData
    
    print "Generating the blank template from the first file"
    # generate blank template Exposure with these characterstics relative to the first exposure:
    # twice as many pixels
    # same RA/Dec center
    # half the scale (degrees/pixel)
    # north up, east to the right
    firstExposure = afwImage.ExposureF()
    firstExposure.readFits(imageDataList[0].filepath)
    templateExposure = makeBlankTemplateExposure(firstExposure)

    print "Building template"
    coaddKaiser.makeBlurredCoadd(templateExposure, imageDataList)
    templateExposure.writeFits(outname)
        
