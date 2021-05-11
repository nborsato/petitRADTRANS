"""
Setup file for package `petitRADTRANS`.
"""
from setuptools import find_packages
from numpy.distutils.core import Extension, setup
import os
import warnings

use_compiler_flags = True

if use_compiler_flags:
    extra_compile_args = ["-O3",
                        "-funroll-loops",
                        "-ftree-vectorize",
                        "-msse",
                        "-msse2",
                        "-m3dnow"]
else:
    extra_compile_args = None


fort_spec = Extension(
    name='petitRADTRANS.fort_spec',
    sources=['petitRADTRANS/fort_spec.f90'],
    extra_compile_args=extra_compile_args)

fort_input = Extension(
    name='petitRADTRANS.fort_input',
    sources=['petitRADTRANS/fort_input.f90'], \
    extra_compile_args=extra_compile_args)

fort_rebin = Extension(
    name='petitRADTRANS.fort_rebin',
    sources=['petitRADTRANS/fort_rebin.f90'], \
    extra_compile_args=extra_compile_args)


extensions = [fort_spec, fort_input, fort_rebin]

def setup_function(extensions):
    setup(name='petitRADTRANS',
          version="2.1.0",
          description='Exoplanet spectral synthesis tool for retrievals',
          long_description=open(os.path.join(
              os.path.dirname(__file__), 'README.rst')).read(),
          long_description_content_tpye='test/x-rst',
          url='https://gitlab.com/mauricemolli/petitRADTRANS',
          author='Paul Mollière',
          author_email='molliere@mpia.de',
          license='MIT License',
          packages=find_packages(),
          include_package_data=True,
          install_requires=['scipy', 'numpy', 'matplotlib', 'h5py'],
          zip_safe=False,
          ext_modules=extensions,
          )

setup_function(extensions)
