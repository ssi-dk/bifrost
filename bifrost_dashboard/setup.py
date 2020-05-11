from setuptools import setup, find_packages

setup(
    name='bifrost-dashboard',
    version='0.0.1',
    description='Dashboard for displaying bifrost information',
    url='https://github.com/ssi-dk/bifrost',
    author="Martin Basterrechea",
    author_email="mbas@ssi.dk",
    packages=find_packages(),
    install_requires=['pymongo'],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'run=bifrost_dashboard.reporter:main'
        ]
    }
)