""" run this from the casadi/debian directory"""

if __name__=='__main__':
    version = '2.0'

    casadi_short = 'numerical optimization / algorithmic differentiation framework'
    casadi_long = '''\
 CasADi is a numerical optimization / algorithmic differentiation framework.
 It can be used from C++ or Python.
 .
 It provides users with a set of building blocks that simplify the process of
 implementing highly efficient gradient-based solvers for numerical optimization
 problems in general and simulation-based nonlinear programs (optimal control)
 in particular. This can be done with a range of different methods including
 direct collocation, direct multiple shooting and indirect methods.
 .
 Contained in the package is a symbolic framework that allows constructing
 large computational graphs made up by matrix-valued operations and a
 state-of-the-art framework for algorithmic differentiation (AD) - also known
 as automatic differentiation - that operates on these computational graphs.
 Also contained is a set of solvers and interfaces to third-party solvers for
 nonlinear programming (NLP), quadratic programming (QP) and initial-value
 problems in ordinary differential equations (ODE) or diffential-algebraic
 equations (DAE). Other features of interest include generation of
 self-contained C-code and symbolic import of models from the Modelica physical
 modeling language and the AMPL algebraic modeling language.
 .
 More information about CasADi can be found on the project website,
 http://casadi.org.
'''
    
    # all the modules we current generate
    stuff = {}
    stuff['core'] = {'dir': 'core'}
    stuff['control'] = {'dir': 'control'}
    #stuff['optimal-control'] = {'dir': 'optimal_control'}
    #stuff['linearsolver-csparse'] = {'dir': 'interfaces/csparse'}
    #stuff['linearsolver-lapack'] = {'dir': 'interfaces/lapack'}
    #stuff['qpsolver-qpoases'] = {'dir': 'interfaces/qpoases'}
    stuff['nlpsolver-snopt'] = {'dir': 'interfaces/snopt'}
    stuff['nlpsolver-ipopt'] = {'dir': 'interfaces/ipopt'}
    #stuff['sdpsolver-dsdp'] = {'dir': 'interfaces/dsdp'}
    #stuff['integrator-sundials'] = {'dir': 'interfaces/sundials'}

    # extra module-specific customization
    stuff['core']['desc_short'] = 'numerical optimization and algorithmic differentiation framework'
    stuff['core']['desc_long'] = '''\
 Core module of CasADi, in particular containing the symbolic framework and a
 self-contained implementation of algorithmic differentiation.
 .
''' + casadi_long

    stuff['nlpsolver-ipopt']['libdeps'] = ['coinor-libipopt1']

    # (end of customization, should need to edit no further)

    control_file = '''\
Source: casadi
Priority: optional
Maintainer: Greg Horn <gregmainland@gmail.com>
Build-Depends: debhelper (>= 9),
               cmake,
               pkg-config,
               coinor-libipopt-dev,
               libblas-dev,
               liblapack-dev,
               swig2.0,
               python-dev,
               python-support (>= 0.4.1),
               python-numpy
Standards-Version: 3.9.5
Section: libs
Homepage: http://casadi.org
Vcs-Git: git://github.com/casadi/casadi.git
Vcs-Browser: https://github.com/casadi/casadi

'''

    for name,info in sorted(stuff.items()):
        directory = info['dir']
        underscores = name.replace('-','_')
        
        # library dependencies
        libdeps = []
        if name != 'core':
            libdeps.append('libcasadi-core'+version+' (= ${binary:Version})')
        if 'libdeps' in info:
            libdeps.extend(info['libdeps'])
        if libdeps:
            libdeps = ',\n         '+',\n         '.join(libdeps)
        else:
            libdeps = ''

        # development library dependencies
        devlibdeps = ['libcasadi-'+name+version+' (= ${binary:Version})']
        if name != 'core':
            devlibdeps = ['libcasadi-core-dev (= ${binary:Version})']
        if 'devlibdeps' in info:
            devlibdeps.extend(info['devlibdeps'])
        if devlibdeps:
            devlibdeps = ',\n         '+',\n         '.join(devlibdeps)
        else:
            devlibdeps = ''

        # short description
        if 'desc_short' in info:
            desc_short = info['desc_short']
        else:
            print name + " missing desc_short"
            desc_short = name + ' module for CasADi optimization framework'
        # long description
        if 'desc_long' in info:
            desc_long = info['desc_long']
        else:
            print name + " missing desc_long"
            desc_long = casadi_long

        # development version descriptions
        dev_desc_short = 'development files for CasADi ' + name + ' module'
        dev_desc_long = desc_long

        ##################### DON'T CHANGE ANYTHING AFTER THIS LINE #####################
        assert len(desc_short) <= 80, (name+ ": short description must be 80 chars or less")
        assert len(dev_desc_short) <= 80, (name+ ": development short description must be 80 chars or less")
        for line in desc_long:
            assert len(line) <= 80, (name+ ": each line of long description must be 80 chars or less")
        for line in dev_desc_long:
            assert len(line) <= 80, (name+ ": each line of development long description must be 80 chars or less")
        control_file += '''\
Package: libcasadi-%(name)s%(version)s
Section: libs
Architecture: any
Depends: ${shlibs:Depends},
         ${misc:Depends}%(libdeps)s
Description: %(desc_short)s
%(desc_long)s

Package: libcasadi-%(name)s-dev
Section: libdevel
Architecture: any
Depends: ${misc:Depends}%(devlibdeps)s
Description: %(dev_desc_short)s
%(dev_desc_long)s

''' % {'name':name,
       'version':version,
       'libdeps':libdeps,
       'devlibdeps':devlibdeps,
       'desc_short':desc_short,
       'dev_desc_short':dev_desc_short,
       'desc_long':desc_long.rstrip(),
       'dev_desc_long':dev_desc_long.rstrip()}

        devfile = '''\
usr/lib/libcasadi_%(name)s.so
usr/lib/pkgconfig/casadi_%(name)s.pc
usr/include/casadi/%(directory)s/*
''' % {'name':underscores,'directory':directory}

        libfile = '''\
usr/lib/libcasadi_%(name)s.so.*
''' % {'name':underscores}

        f = open('libcasadi-'+name+version+'.install', 'w')
        f.write(libfile)

        f = open('libcasadi-'+name+'-dev.install', 'w')
        f.write(devfile)


    control_file += ('''\
Package: python-casadi
Section: python
Architecture: any
Depends: ${python:Depends}, ${shlibs:Depends}, ${misc:Depends},
         libcasadi-core2.0 (= ${binary:Version}),
         python-numpy-abi9,
	 python-numpy (>= 1:1.6.1)
Provides: ${python:Provides}
Description: Python bindings for CasADi
 This package contains Python bindings for CasADi.
 .
''' + casadi_long).rstrip()

    f = open('control', 'w')
    f.write(control_file)
