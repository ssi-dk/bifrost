from setuptools import setup, find_packages

setup(
    name='bifrostlib',
    version='2.0.11',
    description='Datahandling functions for bifrost (later to be API interface)',
    url='https://github.com/ssi-dk/bifrost/tree/master/lib/bifrostlib',
    author="Kim Ng, Martin Basterrechea",
    author_email="kimn@ssi.dk",
    packages=find_packages(),
    install_requires=['ruamel.yaml', 'pymongo', 'pandas', 'numpy', 'python-magic', 'dnspython']
    )
