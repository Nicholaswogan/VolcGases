
from numpy.distutils.core import setup, Extension

sources = ['src/VolcGasesFort.f90',
           'src/minpack/dpmpar.f',
           'src/minpack/enorm.f',
           'src/minpack/fdjac2.f',
           'src/minpack/lmder.f',
           'src/minpack/lmder1.f',
           'src/minpack/lmdif.f',
           'src/minpack/lmpar.f',
           'src/minpack/qrfac.f',
           'src/minpack/qrsolv.f']

extensions = [
Extension(name="VolcGasesFort",
          sources=sources,
          extra_f90_compile_args = ['-O3', '-freal-4-real-8'],
          extra_f77_compile_args = ['-O3', '-freal-4-real-8'],
          f2py_options=['only: solve_gases'])
          ]

setup(name = 'VolcGases',
      packages=['VolcGases'],
      version='2.0',
      ext_modules=extensions)
