from setuptools import setup, find_packages

setup(
    name='fetchtool',
    version='0.0',
    description='A tool to fetch seismograms',
    author='Marcelo Belentani de Bianchi',
    author_email='m.bianchi@iag.usp.br',
    url='https://https://github.com/marcelobianchi/fetchtool',
    packages=find_packages(include=[ 'fetchtool' ]),
    install_requires=[
        'obspy>=1.2.2'
   ],
    extras_require = { },
    setup_requires = [ ],
    tests_require  = [ ],
    entry_points   = { },
    package_data   = { }
)
