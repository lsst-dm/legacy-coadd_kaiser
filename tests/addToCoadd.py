#!/usr/bin/env python
"""
Test lsst.coadd.kaiser.addToCoadd
"""

import os
import math
import pdb # we may want to say pdb.set_trace()
import unittest

import numpy

import eups
import lsst.afw.image as afwImage
import lsst.afw.image.testUtils as imTestUtils
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.pex.exceptions as pexEx
import lsst.coadd.kaiser.kaiserLib as coaddKaiser

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.coadd.kaiser", Verbosity)

dataDir = eups.productDir("afwdata")
if not dataDir:
    raise RuntimeError("Must set up afwdata to run these tests") 

InputMaskedImageName = "871034p_1_MI"
InputMaskedImageNameSmall = "small_MI"

currDir = os.path.abspath(os.path.dirname(__file__))
inFilePath = os.path.join(dataDir, InputMaskedImageName)
inFilePathSmall = os.path.join(dataDir, InputMaskedImageNameSmall)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def referenceAddToCoadd(coadd, depthMap, image, badPixelMask):
    """Reference implementation of lsst.coadd.kaiser.addToCoadd
    
    Unlike lsst.coadd.kaiser.addToCoadd this one does not change the inputs.
    
    Inputs:
    - coadd: coadd before adding image
    - depthMap: depth map before adding image
    - image: masked image to add to coadd
    - badPixelMask: mask of bad pixels to ignore

    Returns two items:
    - coaddArrayList: new coadd as a list of image, mask, variance numpy arrays
    - depthMapArray: new depth map, as a numpy array
    """
    imageArrayList = imTestUtils.arraysFromMaskedImage(image)
    coaddArrayList = imTestUtils.arraysFromMaskedImage(coadd)
    depthMapArray = imTestUtils.arrayFromImage(depthMap)

    badMaskArr = (imageArrayList[1] & badPixelMask) != 0
    for ind in (0, 2):
        coaddArray = coaddArrayList[ind]
        coaddArray += numpy.where(badMaskArr, 0, imageArrayList[ind])
    coaddArray = coaddArrayList[1]
    coaddArray |= numpy.where(badMaskArr, 0, imageArrayList[1])
    depthMapArray += numpy.where(badMaskArr, 0, 1)
    return coaddArrayList, depthMapArray


class addToCoaddTestCase(unittest.TestCase):
    """
    A test case for addToCoadd
    """
    def referenceTest(self, coadd, depthMap, image, badPixelMask):
        """Compare lsst implemenation of addToCoadd to a reference implementation.
        """
        refCoaddArrayList, refDepthMapArray = referenceAddToCoadd(coadd, depthMap, image, badPixelMask)
        coaddKaiser.addToCoadd(coadd, depthMap, image, badPixelMask) # changes the inputs
        coaddArrayList = imTestUtils.arraysFromMaskedImage(coadd)
        maskArr = coaddArrayList[1]
        depthMapArray = imTestUtils.arrayFromImage(depthMap)
        
        for name, ind in (("image", 0), ("mask", 1), ("variance", 2)):
            if not numpy.allclose(coaddArrayList[ind], refCoaddArrayList[ind]):
                errMsgList = (
                    "Computed %s does not match reference for badPixelMask=%s:" % (name, badPixelMask),
                    "computed=  %s" % (coaddArrayList[ind],),
                    "reference= %s" % (refCoaddArrayList[ind],),
                )
                errMsg = "\n".join(errMsgList)
                self.fail(errMsg)
        if numpy.any(depthMapArray != refDepthMapArray):
            errMsgList = (
                "Computed depth map does not match reference for badPixelMask=%s:" % (badPixelMask,),
                "computed=  %s" % (depthMapArray,),
                "reference= %s" % (refDepthMapArray,),
            )
            errMsg = "\n".join(errMsgList)
            self.fail(errMsg)
        
    def testSmall(self):
        """Test addToCoadd on afwdata small image
        """
        image = afwImage.MaskedImageF(inFilePathSmall)
        coadd = afwImage.MaskedImageF(image.getDimensions())
        coadd *= 0.0
        depthMap = afwImage.ImageU(image.getDimensions(), 0)
        for badPixelMask in (0, 0xF, 0xFF, 0xFFF, 0xFFFF):
            self.referenceTest(coadd, depthMap, image, badPixelMask)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(addToCoaddTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
