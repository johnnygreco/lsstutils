#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='lsstutils',
      version='0.1',
      author='Johnny Greco',
      packages=['lsstutils'],
      url='https://github.com/johnnygreco/lsstutils')
