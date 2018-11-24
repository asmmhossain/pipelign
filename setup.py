from setuptools import find_packages, setup
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))

long_description = """
Pipelign
--------

An automated pipeline for multiple sequence alignment.
"""

setup(
    name='pipelign',
    version='0.1',
    author='Mukarram Hossain',
    author_email='asmmh2@cam.ac.uk',
    packages=find_packages(exclude=[dir_path+'/docs']),
    url='https://github.com/asmmhossain/pipelign',
    license='MIT',
    description='A pipeline for automated alignment',
    long_description=long_description,
    install_requires=[
        'biopython',
        'ete3'
    ],
    scripts=[
        dir_path+'/bin/pipelign',
        dir_path+'/bin/gb2fas'
    ]
)
