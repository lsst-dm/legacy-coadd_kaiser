from kaiserLib import *
#from MakeBlurredCoadd import *
# trying to import ANYTHING from MakeBlurredCoadd
# breaks addToMaskedImage unit test with 
# Traceback (most recent call last):
#   File "tests/addToMaskedImage.py", line 67, in testBasics
#     self.referenceTest(blankMaskedImage, inMaskedImage, badPixelMask)
#   File "tests/addToMaskedImage.py", line 43, in referenceTest
#     origOutArrays = imTestUtils.arraysFromMaskedImage(outMaskedImage)
#   File "/Users/rowen/LSST/code/afw-trunk/python/lsst/afw/image/testUtils.py", line 41, in arraysFromMaskedImage
#     arrayFromMask(maskedImage.getMask().get()),
#   File "/Users/rowen/LSST/code/afw-trunk/python/lsst/afw/image/testUtils.py", line 31, in arrayFromMask
#     arr[col, row] = im.getPtr(col, row)
#   File "/Users/rowen/LSST/code/ip_diffim-trunk/python/lsst/ip/diffim/diffimLib.py", line 1991, in <lambda>
#     __getattr__ = lambda self, name: _swig_getattr(self, MaskU, name)
#   File "/Users/rowen/LSST/code/ip_diffim-trunk/python/lsst/ip/diffim/diffimLib.py", line 40, in _swig_getattr
#     raise AttributeError,name
# AttributeError: getPtr
