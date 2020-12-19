'''Setup script for ATACFragQC'''

from __future__ import absolute_import

import os.path
from setuptools import setup
from ATACFragQC import __version__

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, 'README.md')) as fid:
    README = fid.read()

# This call to setup() does all the work
setup(
    name='ATACFragQC',
    version=__version__,
    description='Fragment Quality Control for ATAC-seq',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/0CBH0/ATACFragQC',
    author='0CBH0',
    author_email='maodatou88@163.com',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    packages=['ATACFragQC'],
    include_package_data=True,
    install_requires=[
        'pysam', 'pandas', 'numpy', 'plotnine', 'Pillow'
    ],
    python_requires='>=3.6',
    entry_points={'console_scripts': ['ATACFragQC=ATACFragQC.__main__:main']},
)
