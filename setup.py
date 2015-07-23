#from distutils.core import setup
#from distutils.extension import Extension
#from distutils.command.sdist import sdist as _sdist

from setuptools import setup
from setuptools import Extension
from setuptools.command.sdist import sdist as _sdist
from fmrc import util

#borrowed from online code: http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
try:
    from Cython.Distutils import build_ext
except ImportError:
    useCython = False
else:
    useCython = True

import numpy as np

cmdClass = {}
extModules = []

if useCython:
    extModules += [Extension('fmrc.correct', ['fmrc/correct.pyx'])]
    cmdClass.update({'build_ext': build_ext})

    #this is also from the stackoverflow link above, used to auto-compile when you do the sdist command
    class sdist(_sdist):
        def run(self):
            # Make sure the compiled Cython files in the distribution are up-to-date
            from Cython.Build import cythonize
            cythonize('fmrc/correct.pyx', include_path=[np.get_include()])
            _sdist.run(self)
    cmdClass['sdist'] = sdist

else:
    extModules += [Extension('fmrc.correct', ['fmrc/correct.c'])]

setup(name='fmrc',
      version=util.VERSION,
      description='Corrects errors in short reads from high-throughput sequencing',
      url='http://github.com/sgreenstein/fmrc',
      author='Seth Greenstein',
      author_email='sgreens@cs.unc.edu',
      license='MIT',
      install_requires=['pysam', 'numpy'],
      scripts=['bin/fmrc'],
      packages=['fmrc'],
      zip_safe=False,
      include_dirs=[np.get_include()],
      ext_modules=extModules,
      cmdclass=cmdClass)
