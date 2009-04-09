#!/usr/bin/env python
"""
Test lsst.coadd.kaiser.medianBinapprox

@todo: expand tests to lots of other data including some unsual distributions
"""

import os
import math
import pdb # we may want to say pdb.set_trace()
import unittest

import numpy

import eups
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexEx
import lsst.coadd.kaiser as coaddKaiser
import lsst.afw.image.testUtils as imTestUtils

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.coadd.kaiser", Verbosity)

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests") 

InputImageNameSmall = "small_MI_img.fits"

currDir = os.path.abspath(os.path.dirname(__file__))
inFilePathSmallImage = os.path.join(dataDir, InputImageNameSmall)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def refMedian(inArr):
    """Compute the median of an array of any shape.
    """
    nPts = inArr.size
    if nPts == 0:
        raise RuntimeError("image array has no pixels")
    linArr = inArr.flatten()
    linArr.sort()
    ctrInd = (nPts-1)//2
    if nPts % 2 == 0:
        # array has even length; return mean of two central values
        return (linArr[ctrInd] + linArr[ctrInd+1]) / 2.0
    else:
        # array has odd length; return central value
        return linArr[ctrInd]
        

class medianBinApproxTestCase(unittest.TestCase):
    def _testOneImage(self, image, nBins):
        """Test medianBinapproxImage on one image
        """
        med = coaddKaiser.medianBinapproxImage(image, int(nBins))
        imArr = imTestUtils.arrayFromImage(image)
        refMed = refMedian(imArr)
        stdDev = imArr.std()
        err = med - refMed
        maxErr = stdDev / float(nBins)
#        print "med=%s, refMed=%s, err=%s, maxErr=%s" % (med, refMed, err, maxErr)
        if abs(err) > maxErr:
            self.fail("Computed median error too large")

    """
    A test case for the medianBinApprox Class
    """
    def testSmallImage(self):
        """Test median of small image
        """
        image = afwImage.ImageF(inFilePathSmallImage)
        for nBins in (10, 100, 1000):
            self._testOneImage(image, nBins)
    
        

         
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(medianBinApproxTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
