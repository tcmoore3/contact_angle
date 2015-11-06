from setuptools import setup
from setuptools.command.test import test as TestCommand
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running contact_angle requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)

requirements = [line.strip() for line in open('requirements.txt').readlines()]


class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(['contact_angle'])
        sys.exit(errcode)


setup(name='contact_angle',
      version='beta',
      description='Routines for calculating contact angles',
      url='http://github.com/tcmoore3/contact_angle',
      author='Timothy C. Moore',
      author_email='timothy.c.moore@vanderbilt.edu',
      license='MIT',
      packages=['contact_angle'],
      install_requires=requirements,
      zip_safe=False,
      test_suite='tests',
      cmdclass={'test': PyTest},
      extras_require={'utils': ['pytest']},
)
