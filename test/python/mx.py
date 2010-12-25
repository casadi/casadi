from casadi import *
from numpy import *
import unittest
from types import *

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
       	pass

    def test_MXfunction1(self):
	# check if x->2*x
	# evaluates correctly for x=3
	x = MX("x")
	y = 2*x
	f = MXFunction([x],[y])
	self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
	self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")

        f.init()
	f.setInput(3,0);
	f.evaluate()
	yt = f.getOutput()
	self.assertEqual(type(yt),TupleType,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(len(yt),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	y=yt[0]
	self.assertEqual(type(y),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertAlmostEqual(y, 2*3,10)

    def test_MXfunction2(self):
	# check if [x,y]->[y+x,y*x]
	# evaluates correctly for x=3,y=7
	x = MX("x")
	y = MX("y")
	f = MXFunction([x,y],[x+y,y*x])
	self.assertEqual(f.getNumInputs(),2,"MXFunction fails to indicate correct number of inputs")
	self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")

        f.init()
	f.setInput(3,0);
	f.setInput(7,1);
	f.evaluate()
	zt1 = f.getOutput(0)
	zt2 = f.getOutput(1)
	self.assertEqual(type(zt1),TupleType,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(type(zt2),TupleType,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(len(zt1),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	self.assertEqual(len(zt2),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	z1=zt1[0]
	z2=zt2[0]
	self.assertEqual(type(z1),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(type(z2),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertAlmostEqual(z2, 21,10)
	self.assertAlmostEqual(z1, 10,10)



    def test_MXfunction3(self):
	# check if [x,y]->[y+x,y*x]
	# evaluates correctly for x=3,y=7
	# now with single input, multi output
	xy = MX("xy",2)
	f = MXFunction([xy],[xy[0]+xy[1],xy[0]*xy[1]])
	self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
	self.assertEqual(f.getNumOutputs(),2,"MXFunction fails to indicate correct number of outputs")
        f.init()
	f.setInput([3,7],0);
	f.evaluate()
	zt1 = f.getOutput(0)
	zt2 = f.getOutput(1)
	self.assertEqual(type(zt1),TupleType,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(type(zt2),TupleType,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(len(zt1),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	self.assertEqual(len(zt2),1,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	z1=zt1[0]
	z2=zt2[0]
	self.assertEqual(type(z1),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(type(z2),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertAlmostEqual(z2, 21,10)
	self.assertAlmostEqual(z1, 10,10)

    def test_MXfunction4(self):
	# check if [x,y]->[y+x,y*x]
	# evaluates correctly for x=3,y=7
	# now with single input, single output
	xy = MX("xy",2)
	z=MX("z",2)
	z[0]=xy[0]+xy[1]
	z[1]=xy[0]*xy[1]
	f = MXFunction([xy],[z])
	self.assertEqual(f.getNumInputs(),1,"MXFunction fails to indicate correct number of inputs")
	self.assertEqual(f.getNumOutputs(),1,"MXFunction fails to indicate correct number of outputs")
        f.init()
	f.setInput([3,7],0);
	f.evaluate()
	zt = f.getOutput()
	self.assertEqual(type(zt),TupleType,"Output of MXFunction is expected to be tuple of floats")
	z1=zt[0]
	z2=zt[1]
	self.assertEqual(len(zt),2,"Output of MXFunction was tuple of floats, as expected, but length is incorrect.")
	self.assertEqual(type(z1),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertEqual(type(z2),float,"Output of MXFunction is expected to be tuple of floats")
	self.assertAlmostEqual(z2, 21,10)
	self.assertAlmostEqual(z1, 10,10)

if __name__ == '__main__':
    unittest.main()

