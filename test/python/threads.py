#
#     This file is part of CasADi.
#
#     CasADi -- A symbolic framework for dynamic optimization.
#     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
#                             KU Leuven. All rights reserved.
#     Copyright (C) 2011-2014 Greg Horn
#
#     CasADi is free software; you can redistribute it and/or
#     modify it under the terms of the GNU Lesser General Public
#     License as published by the Free Software Foundation; either
#     version 3 of the License, or (at your option) any later version.
#
#     CasADi is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     Lesser General Public License for more details.
#
#     You should have received a copy of the GNU Lesser General Public
#     License along with CasADi; if not, write to the Free Software
#     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
#
from casadi import *
import casadi as ca
import numpy
import unittest
from types import *
from helpers import *
import pickle
import os
import sys
import datetime
import threading
import time
import multiprocessing

target_print_thread = 0.25

class PrintNowThread:
    def __init__(self):
        self._thread = threading.Thread(target=self.run)
        self._thread.daemon = True
        self.timestamps = []

    def start(self):
        self._thread.start()

    def run(self):
        while True:
            print("____________________________________________________________________")
            print(
                "                                                      The time is",
                datetime.datetime.now(),
            )
            print("____________________________________________________________________")
            time.sleep(target_print_thread)
            self.timestamps.append(time.time())
            
class Threadstests(casadiTestCase):

  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_GIL_release_wall_time(self):




        timerthread = PrintNowThread()
        timerthread.start()


        print("foo")
        import casadi as ca
        import numpy as np

        N = 400  # number of control intervals

        opti = ca.Opti()  # Optimization problem

        # ---- decision variables ---------
        X = opti.variable(2, N + 1)  # state trajectory
        pos = X[0, :]
        speed = X[1, :]
        U = opti.variable(1, N)  # control trajectory (throttle)
        T = opti.variable()  # final time

        # ---- objective          ---------
        opti.minimize(T)  # race in minimal time

        # ---- dynamic constraints --------
        f = lambda x, u: ca.vertcat(x[1], u - x[1])  # dx/dt = f(x,u)

        dt = T / N  # length of a control interval
        for k in range(N):  # loop over control intervals
            # Runge-Kutta 4 integration
            k1 = f(X[:, k], U[:, k])
            k2 = f(X[:, k] + dt / 2 * k1, U[:, k])
            k3 = f(X[:, k] + dt / 2 * k2, U[:, k])
            k4 = f(X[:, k] + dt * k3, U[:, k])
            x_next = X[:, k] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            opti.subject_to(X[:, k + 1] == x_next)  # close the gaps

        # ---- path constraints -----------
        limit = lambda pos: 1 - ca.sin(2 * ca.pi * pos) / 2
        opti.subject_to(speed <= limit(pos))  # track speed limit
        opti.subject_to(opti.bounded(0, U, 1))  # control is limited

        # ---- boundary conditions --------
        opti.subject_to(pos[0] == 0)  # start at position 0 ...
        opti.subject_to(speed[0] == 0)  # ... from stand-still
        opti.subject_to(pos[-1] == 1)  # finish line at position 1

        # ---- misc. constraints  ----------
        opti.subject_to(T >= 0)  # Time must be positive

        # ---- initial values for solver ---
        opti.set_initial(speed, 1)
        opti.set_initial(T, 1)

        # ---- solve NLP              ------
        opti.solver("ipopt")  # set numerical backend

        fun = opti.to_function('opti',[T],[X])
        raw_solver = fun.find_function("solver")

        print("fun(1)")
        fun(1)
        print("done")

        base_time = fun.stats()["t_wall_nlp_jac_g"]

        print(base_time)

        n_thread = 2

        from multiprocessing.pool import ThreadPool

        buffers = [fun.buffer() for _ in range(n_thread)]


        def casadi_calc(i,**kwargs):
            [buffer,trigger] = buffers[i]
            
            T = np.ones((1,1))
            X = np.random.random((2, N + 1))
            buffer.set_arg(0, memoryview(T))
            buffer.set_res(0, memoryview(X))
            
            trigger()


        p = ThreadPool(n_thread)

        results = p.map(casadi_calc,range(n_thread))

        results = [raw_solver.stats(i+1)["t_wall_nlp_jac_g"] for i in range(n_thread)]
        print(base_time,results)


        # Evaluation wall time should stay the same with two threads
        for e in results:
            print(e/base_time)
            self.assertTrue(0.7 <= e/base_time <= 1.3)

        print(timerthread.timestamps)
        for dt in np.diff(np.array(timerthread.timestamps)):
            print("target_print_thread dt",dt)
            self.assertTrue(target_print_thread*0.9<=dt<=target_print_thread*1.1)
        
  @memory_heavy()
  @requires_nlpsol("ipopt")
  def test_GIL_release_stress_test(self):

    timerthread = PrintNowThread()
    timerthread.start()


    print("foo")
    import casadi as ca
    import numpy as np

    N = 100  # number of control intervals

    opti = ca.Opti()  # Optimization problem

    # ---- decision variables ---------
    X = opti.variable(2, N + 1)  # state trajectory
    pos = X[0, :]
    speed = X[1, :]
    U = opti.variable(1, N)  # control trajectory (throttle)
    T = opti.variable()  # final time


    S = ca.SX.sym("S",4,4)

    callback_buffer = []

    class MyCallback(ca.Callback):
      def __init__(self, name, opts={}):
        ca.Callback.__init__(self)
        self.construct(name, opts)
        self.counter = 0



      def eval(self,argin):
        # Purposefully not thread-safe to solicit crashes
        callback_buffer.append(argin)
        time.sleep(0.001)
        self.counter = self.counter + 1
        #print(ca.det(S))
        if argin[0]>10:
            raise Exception("let's try an exception")
        
        return [argin[0]]
        
    mycallback = MyCallback("mycallback",{"enable_fd":True,"fd_method":"backward"})

    print(mycallback(T))


    # ---- objective          ---------
    opti.minimize(mycallback(T))  # race in minimal time

    # ---- dynamic constraints --------
    f = lambda x, u: ca.vertcat(x[1], u - x[1])  # dx/dt = f(x,u)

    dt = T / N  # length of a control interval
    for k in range(N):  # loop over control intervals
        # Runge-Kutta 4 integration
        k1 = f(X[:, k], U[:, k])
        k2 = f(X[:, k] + dt / 2 * k1, U[:, k])
        k3 = f(X[:, k] + dt / 2 * k2, U[:, k])
        k4 = f(X[:, k] + dt * k3, U[:, k])
        x_next = X[:, k] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        opti.subject_to(X[:, k + 1] == x_next)  # close the gaps

    # ---- path constraints -----------
    limit = lambda pos: 1 - ca.sin(2 * ca.pi * pos) / 2
    opti.subject_to(speed <= limit(pos))  # track speed limit
    opti.subject_to(opti.bounded(0, U, 1))  # control is limited

    # ---- boundary conditions --------
    opti.subject_to(pos[0] == 0)  # start at position 0 ...
    opti.subject_to(speed[0] == 0)  # ... from stand-still
    opti.subject_to(pos[-1] == 1)  # finish line at position 1

    # ---- misc. constraints  ----------
    opti.subject_to(T >= 0)  # Time must be positive

    # ---- initial values for solver ---
    opti.set_initial(speed, 1)
    opti.set_initial(T, 1)


    foo = []

    class IterationCallback(ca.Callback):
      def __init__(self,nx, ng, np):
        ca.Callback.__init__(self)
        self.foo = []

        self.nx = nx
        self.ng = ng
        self.np = np
        self.construct("mycallback", {})

      def get_n_in(self): return nlpsol_n_out()
      def get_n_out(self): return 1


      def get_sparsity_in(self, i):
        n = nlpsol_out(i)
        if n=='f':
          return Sparsity. scalar()
        elif n in ('x', 'lam_x'):
          return Sparsity.dense(self.nx)
        elif n in ('g', 'lam_g'):
          return Sparsity.dense(self.ng)
        else:
          return Sparsity(0,0)
      def eval(self, arg):
        # Purposefully not thread-safe to solicit crashes
        print("yay")
        try:
            raise Exception()
        except:
            print("caught")
        foo.append(arg)
        return [0]
    iteration_callback = IterationCallback(opti.nx, opti.np, opti.ng)
        
    # ---- solve NLP              ------
    opti.solver("ipopt",{"iteration_callback": iteration_callback,"ipopt.hessian_approximation":"limited-memory"})  # set numerical backend


    fun = opti.to_function('opti',[T],[X])
    raw_solver = fun.find_function("solver")

    n_thread = 4

    from multiprocessing.pool import ThreadPool

    buffers = [fun.buffer() for _ in range(n_thread)]


    def casadi_calc(i,**kwargs):
        [buffer,trigger] = buffers[i]
        
        T = np.ones((1,1))
        X = np.random.random((2, N + 1))
        buffer.set_arg(0, memoryview(T))
        buffer.set_res(0, memoryview(X))
        
        trigger()
        
        return X


    p = ThreadPool(n_thread)

    results = p.map(casadi_calc,range(n_thread))
    print(results)
    
    for i in range(n_thread):
        self.checkarray(results[0],results[i],digits=14)
        self.assertTrue(raw_solver.stats(i+1)["success"])

    print(timerthread.timestamps)
    for dt in np.diff(np.array(timerthread.timestamps)):
        print("target_print_thread dt",dt)
        self.assertTrue(target_print_thread*0.9<=dt<=target_print_thread*1.1)
    
    [raw_solver.stats(i+1)["t_wall_nlp_jac_g"] for i in range(n_thread)]

  def test_threadsafe_symbolics(self):
  
    if "CASADI_WITH_THREADSAFE_SYMBOLICS" not in CasadiMeta.compiler_flags(): return
        
    timerthread = PrintNowThread()
    timerthread.start()




    x = ca.SX.sym("x")

    x = 2*x
    x = 2*x

    x = 0

    n_thread = 4

    from multiprocessing.pool import ThreadPool

    A = ca.SX.sym("A",5,5)

    def casadi_calc(i,**kwargs):

        for i in range(10):
            ca.det(A)
            
        # Circular dependency
        foo = {"a": A}
        bar = {"bar":foo}
        foo["baz"] = bar

        N = 50  # number of control intervals

        opti = ca.Opti()  # Optimization problem

        # ---- decision variables ---------
        X = opti.variable(2, N + 1)  # state trajectory
        pos = X[0, :]
        speed = X[1, :]
        U = opti.variable(1, N)  # control trajectory (throttle)
        T = opti.variable()  # final time


        callback_buffer = []

        class MyCallback(ca.Callback):
          def __init__(self, name, opts={}):
            ca.Callback.__init__(self)
            self.construct(name, opts)

          def eval(self,argin):
            # Purposefully not thread-safe to solicit crashes
            print("yay")
            # Purposefully not thread-safe to solicit crashes
            callback_buffer.append(argin)
            time.sleep(0.01)
            return [argin[0]**2]
            
        mycallback = MyCallback("mycallback",{"enable_fd":True})

        print(mycallback(T))


        # ---- objective          ---------
        opti.minimize(mycallback(T))  # race in minimal time

        # ---- dynamic constraints --------
        f = lambda x, u: ca.vertcat(x[1], u - x[1])  # dx/dt = f(x,u)

        dt = T / N  # length of a control interval
        for k in range(N):  # loop over control intervals
            # Runge-Kutta 4 integration
            k1 = f(X[:, k], U[:, k])
            k2 = f(X[:, k] + dt / 2 * k1, U[:, k])
            k3 = f(X[:, k] + dt / 2 * k2, U[:, k])
            k4 = f(X[:, k] + dt * k3, U[:, k])
            x_next = X[:, k] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            opti.subject_to(X[:, k + 1] == x_next)  # close the gaps

        # ---- path constraints -----------
        limit = lambda pos: 1 - ca.sin(2 * ca.pi * pos) / 2
        opti.subject_to(speed <= limit(pos))  # track speed limit
        opti.subject_to(opti.bounded(0, U, 1))  # control is limited

        # ---- boundary conditions --------
        opti.subject_to(pos[0] == 0)  # start at position 0 ...
        opti.subject_to(speed[0] == 0)  # ... from stand-still
        opti.subject_to(pos[-1] == 1)  # finish line at position 1

        # ---- misc. constraints  ----------
        opti.subject_to(T >= 0)  # Time must be positive

        # ---- initial values for solver ---
        opti.set_initial(speed, 1)
        opti.set_initial(T, 1)

        # ---- solve NLP              ------
        opti.solver("ipopt")  # set numerical backend

        fun = opti.to_function('opti',[T],[X])
        
        fun(1)
        
        return True
    
    


    p = ThreadPool(n_thread)

    results = p.map(casadi_calc,range(n_thread))
  
if __name__ == '__main__':
    unittest.main()
