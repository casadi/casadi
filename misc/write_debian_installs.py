""" run this from the casadi/debian directory"""

if __name__=='__main__':
    version = '2.0'

    # all the modules we current generate
    stuff = {}
    stuff['core'] = {'dir': 'core'}
    stuff['control'] = {'dir': 'control'}
    #stuff['convex-programming'] = {'dir': 'convex_programming'}
    #stuff['csparse-interface'] = {'dir': 'interfaces/csparse_interface'}
    #stuff['dsdp-interface'] = {'dir': 'interfaces/dsdp_interface'}
    stuff['integration'] = {'dir': 'integration'}
    stuff['ipopt-interface'] = {'dir': 'interfaces/ipopt'}
    #stuff['lapack-interface'] = {'dir': 'interfaces/lapack'}
    stuff['nonlinear-programming'] = {'dir': 'nonlinear_programming'}
    #stuff['optimal-control'] = {'dir': 'optimal_control'}
    #stuff['qpoases-interface'] = {'dir': 'interfaces/qpoases'}
    stuff['snopt-interface'] = {'dir': 'interfaces/snopt'}
    #stuff['sundials-interface'] = {'dir': 'interfaces/sundials'}

    # extra module-specific customization
    stuff['core']['desc_short'] = 'CasADi core module'
    stuff['core']['desc_long'] = '''\
 Core module of CasADi, in particular containing the symbolic framework and a self-contained implementation
 of algorithmic differentiation.
 .
 CasADi is a symbolic framework for algorithmic differentiation and numerical optimization.
 For more information, see http://casadi.org.
'''

    stuff['ipopt-interface']['libdeps'] = ['coinor-libipopt1']

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
               swig,
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
            desc_short = 'CasADi module: ' + name + '. No further information available.'

        # long description
        if 'desc_long' in info:
            desc_long = info['desc_long']
        else:
            print name + " missing desc_long"
            desc_long = '''\
 CasADi is a symbolic framework for algorithmic differentiation and numerical optimization.
 For more information, see http://casadi.org.
'''
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
Description: development files for CasADi %(name)s
 .
 CasADi is a symbolic framework for algorithmic differentiation and numerical optimization.
 For more information, see http://casadi.org.
''' % {'name':name,
       'version':version,
       'libdeps':libdeps,
       'devlibdeps':devlibdeps,
       'desc_short':desc_short,
       'desc_long':desc_long.rstrip()}

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


    control_file += '''\
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
 CasADi is a symbolic framework for algorithmic differentiation and numerical optimization.
 For more information, see http://casadi.org.
'''.rstrip()

    f = open('control', 'w')
    f.write(control_file)
