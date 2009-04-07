"""
Test lsst.coadd.kaiser.CoaddComponent
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

class CoaddComponentTestCase(unittest.TestCase):
    """
    A test case for the CoaddComponent Class
    """
    def testBasics(self):
        """
        Make sure swig interface works
        """
        testExposure = afwImage.ExposureF(inFilePathSmall)
        blankImage = afwImage.ImageD(5, 5)
        blankImage.set(0)
        nullKernel = afwMath.DeltaFunctionKernel(2, 2, 5, 5)
        print "testExposure=%r" % (testExposure,)
        coaddComp = coaddKaiser.CoaddComponent(testExposure, nullKernel)
        #... now verify that we got the right component back

         
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CoaddComponentTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
