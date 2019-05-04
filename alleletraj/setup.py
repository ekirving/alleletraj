#!/usr/bin/env python
# -*- coding: utf-8 -*-

# standard modules
from setuptools import setup, find_packages

setup(name='alleletraj',
      version='1.0',
      description='Selection trajectories of genetic variants underlying domestic animal traits',
      url='https://github.com/ekirving/alleletraj',
      author='Evan K. Irving-Pease',
      author_email='evan.irvingpease@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
            # TODO add these
      ],
      include_package_data=True,
      zip_safe=False)
