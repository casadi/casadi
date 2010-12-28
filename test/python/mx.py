from casadi import *
from numpy import *
import unittest
from types import *

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
       	pass

    def test_MX1(self):
	x = MX("x",2,3)
	self.assertEqual(x.size1(),2,"MX fails to indicate its size1")
	self.assertEqual(x.size2(),3,"MX fails to indicate its size2")

    def test_MXvertcat(self):
	x = MX("x",1,3)
	y = MX("y",1,3)
	z=vertcat(x,y)
	self.assertEqual(z.size1(),2,"MX fails to indicate its size1")
	self.assertEqual(z.size2(),3,"MX fails to indicate its size2")

    def test_MXindices(self):
	x = MX("x",2,3)
	for i in range(2):
		for j in range(3):
			a = "x(%d,%d)" % (i,j)
			b = str(x[i,j])
			self.assertEqual(a,b,"MX indexing is mixed up")

    def test_MXslicing(self):
	x = MX("x",2,3)
	y=x.getColumn(0)
	self.assertEqual(y.size1(),2,"GetColumn returns MX of wrong dimension")
	self.assertEqual(y.size2(),1,"GetColumn returns MX of wrong dimension")

	self.assertEqual(str(y[0]),"[[x(0,0),];[x(1,0),];][0]","Slicing error")
	self.assertEqual(str(y[1]),"[[x(0,0),];[x(1,0),];][1]","Slicing error")

	z=x.getRow(1)
	self.assertEqual(z.size1(),1,"GetRow returns MX of wrong dimension")
	self.assertEqual(z.size2(),3,"GetRow returns MX of wrong dimension")

	self.assertEqual(str(z[0]),"[[x(1,0),x(1,1),x(1,2),];](0,0)","Slicing error")
	self.assertEqual(str(z[1]),"[[x(1,0),x(1,1),x(1,2),];](0,1)","Slicing error")
	self.assertEqual(str(z[2]),"[[x(1,0),x(1,1),x(1,2),];](0,2)","Slicing error")

	x.getRow(0)
	x.getColumn(0)
	x.getColumn(1)
	x.getColumn(2)

    def test_MX2(self):
	U = MX("U",10,2)
	u = U.getRow(1)

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

    def test_MXfunction3b(self):
	# check if [x,y]->[y+x,y*x]
	# evaluates correctly for x=3,y=7
	# now with single input, multi output
	xy = MX("xy",1,2)
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
	z=vertcat([xy[0]+xy[1],xy[0]*xy[1]])
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

    def test_MXfunction5(self):
	# check if [x,y]->[y+x,y*x]
	# evaluates correctly for x=3,y=7
	# now with single input, single output
	xy = MX("xy",2)
	z=horzcat([xy[0]+xy[1],xy[0]*xy[1]])
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

