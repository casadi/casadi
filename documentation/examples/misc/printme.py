from casadi import *



#! Note: the actual printing doesn't show up in this pdf, due to limitations of the reporting software.
#! Better run this example in a terminal.



a = SX("a")
b = SX("b")

c = a+b
c = c.printme(13)

d = c**2

print d

f = SXFunction([a,b],[d])
f.init()
f.input(0).set(4)
f.input(1).set(3)

#! When the graph is evaluated, a printout of c will occur (if you have set WITH_PRINTME to ON in CMakeCache.txt)
#! Printout reads '|> 13: 7'
#! 13 is an identifier of choice, 7 is the numerical value of c
f.evaluate()

J = f.jacobian(0,0)
J.init()

J.init()
J.input(0).set(2)
J.input(1).set(9)

#! The first derivative still depends on c
#! Printout reads '|> 13: 11'
J.evaluate()


J = J.jacobian(0,0)
J.init()

J.init()
J.input(0).set(2)
J.input(1).set(9)
#! second derivative doesn't, so we don't get a printout
J.evaluate()
