from setuptools import setup, find_packages

VERSION = '1.3.3'

setup(
    name='cytopus',
    version=VERSION,
    author="Thomas Walle",
    packages=find_packages(),
    install_requires = [
        #"pandas>1.3",
        #"numpy>1.2",
        "networkx>2.7",
        #"matplotlib>3.4"
        ],
    include_package_data=True,
    package_data={'cytopus': ['data/*.txt']},
)
