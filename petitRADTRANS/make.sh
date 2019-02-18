f2py -c --opt='-O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow' -m fort_input fort_input.f90
f2py -c --opt='-O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow' -m fort_spec fort_spec.f90
f2py -c --opt='-O3 -funroll-loops -ftree-vectorize -ftree-loop-optimize -msse -msse2 -m3dnow' -m fort_rebin fort_rebin.f90
