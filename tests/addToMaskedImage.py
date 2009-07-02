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

class addToMaskedImageTestCase(unittest.TestCase):
    """
    A test case for addToCoadd
    """
    def referenceTest(self, outMaskedImage, inMaskedImage, badPixelMask):
        """Compare lsst implemenation of addToCoadd to a reference implementation.
        """
        origOutArrays = imTestUtils.arraysFromMaskedImage(outMaskedImage)
        inArrays = imTestUtils.arraysFromMaskedImage(inMaskedImage)
        coaddKaiser.addToCoadd(outMaskedImage, inMaskedImage, badPixelMask)
        computedOutArrays = imTestUtils.arraysFromMaskedImage(outMaskedImage)
        
        badMaskArr = (inArrays[1] & badPixelMask) != 0
        refOutArrays = (
            numpy.where(badMaskArr, origOutArrays[0], origOutArrays[0] + inArrays[0]),
            numpy.where(badMaskArr, origOutArrays[1], origOutArrays[1] | inArrays[1]),
            numpy.where(badMaskArr, origOutArrays[2], origOutArrays[2] + inArrays[2]),
        )
        for name, ind in (("image", 0), ("mask", 1), ("variance", 2)):
            if not numpy.allclose(computedOutArrays[ind], refOutArrays[ind]):
                errMsgList = (
                    "Computed %s does not match reference for badPixelMask=%s:" % (name, badPixelMask),
                    "input=     %s" % (inArrays[ind],),
                    "orig out=  %s" % (origOutArrays[ind],),
                    "computed=  %s" % (computedOutArrays[ind],),
                    "reference= %s" % (refOutArrays[ind],),
                    "in mask=   %s" % (inArrays[2],),
                    "badMaskArr=%s" % (badMaskArr,),
                )
                errMsg = "\n".join(errMsgList)
                self.fail(errMsg)
        
        
        
    def testSmall(self):
        """Test addToCoadd on afwdata small image
        """
        inMaskedImage = afwImage.MaskedImageF(inFilePathSmall)
        outMaskedImage = afwImage.MaskedImageF(inMaskedImage.getWidth(), inMaskedImage.getHeight())
        for badPixelMask in (0, 0xF, 0xFF, 0xFFF, 0xFFFF):
            self.referenceTest(outMaskedImage, inMaskedImage, badPixelMask)

         
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(addToMaskedImageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
