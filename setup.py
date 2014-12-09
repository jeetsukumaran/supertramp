#! /usr/bin/env python

#from distutils.core import setup
from setuptools import setup

setup(
    name="supertramp",
    version="0.1.0",
    author="Jeet Sukumaran",
    author_email="jeetsukumaran@gmail.com",
    packages=["supertramp", "test"],
    scripts=["bin/supertramp-simulate.py",
            "bin/supertramp-generate-jobs.py",
            "bin/supertramp-summarize.py",
            "bin/supertramp-summarize-trees.py",
            ],
    url="http://pypi.python.org/pypi/supertramp/",
    license="LICENSE.txt",
    description="A Project",
    long_description=open("README.rst").read(),
    # install_requires=[ ],
)
