gfortran src/VolcGasesFort.f90 src/minpack/dpmpar.f src/minpack/enorm.f \
src/minpack/fdjac2.f src/minpack/lmder.f src/minpack/lmder1.f \
src/minpack/lmdif.f src/minpack/lmpar.f src/minpack/qrfac.f \
src/minpack/qrsolv.f src/VolcGasesFort_wrapper.f90 -shared -fPIC -o VolcGasesFort.so -O3

rm volc.mod volc_wrapper.mod