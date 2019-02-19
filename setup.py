

from distutils.core import setup

setup(
    name='SPECIES-XPS',
    version='0.1.0',
    author='N. Johansson',
    author_email='niclas_v.johansson@maxiv.lu.se',
    packages=['species-xps'],
    url='http://github.com/Smisklas/species-xps/',
    license='LICENSE.txt',
    description='Library for the processing of XPS data',
    long_description=open('README.txt').read(),
    install_requires=[
        "datetime",
        "numpy",
        're',
    ],
)
