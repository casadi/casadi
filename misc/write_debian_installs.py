""" run this from the casadi/debian directory"""

if __name__=='__main__':
    version = '2.0'

    for (name,directory) in [('core','core'),
                             ('control','control'),
                             ('optimal-control','optimal_control'),
                             ('convex-programming','convex_programming'),
                             ('integration','integration'),
                             ('nonlinear-programming','nonlinear_programming'),
                             ('ipopt-interface','interfaces/ipopt'),
                             ('sundials-interface','interfaces/sundials')]:
        underscores = name.replace('-','_')
        libfile = '''\
usr/lib/libcasadi_%(name)s.so.*
''' % {'name':underscores}

        devfile = '''\
usr/lib/libcasadi_%(name)s.so
usr/lib/pkgconfig/casadi_%(name)s.pc
usr/include/casadi/%(directory)s/*
''' % {'name':underscores,'directory':directory}

        f = open('libcasadi-'+name+version+'.install', 'w')
        f.write(libfile)

        f = open('libcasadi-'+name+'-dev.install', 'w')
        f.write(devfile)
