from setuptools import setup, find_packages

setup(
    name='bifrost-dashboard',
    version='0.0.6',
    description='Dashboard for displaying bifrost information',
    url='https://github.com/ssi-dk/bifrost',
    author="Martin Basterrechea",
    author_email="mbas@ssi.dk",
    packages=find_packages(),
    package_data={
        "bifrost_dashboard": ['data/assets/*.css', 'data/assets/*.js',
                              'data/assets/img/*.png', 'data/assets/img/*.svg']
    },
    include_package_data=True,
    install_requires=[
        # 'pymongo',
        'pandas',
        'dash',
        'bifrostapi==0.1.1',
        'Flask-Caching',
        'dash-auth',
        'requests',
        'dash-bootstrap-components',
        'pyyaml',
        #'uwsgi'  # doesnt work for me
        ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'bifrost_dashboard_dev=bifrost_dashboard.reporter:main'
        ]
    },
    scripts=['bin/bifrost_dashboard_start']
)
