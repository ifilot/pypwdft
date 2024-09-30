from setuptools import Extension, setup
import os
import sys
import re

PKG = "pypwdft"
VERSIONFILE = os.path.join(os.path.dirname(__file__), PKG, "_version.py")
verstr = "unknown"
try:
    verstrline = open(VERSIONFILE, "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        print(r"Unable to find version in %s" % (VERSIONFILE,))
        raise RuntimeError(r"If %s.py exists, it is required to be well-formed" % (VERSIONFILE,))

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

setup(
    name=PKG,
    version=verstr,
    author="Ivo Filot",
    author_email="ivo@ivofilot.nl",
    description="Python package for performing plane-wave DFT calculations",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://www.github.com/ifilot/pypwdft",
    packages=[PKG],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
    ],
    python_requires='>=3.5',
    zip_safe=False
)
