from setuptools import setup, find_packages

import geoplot


def get_content(filename):
    with open(filename, "r") as fh:
        return fh.read()


requirements = get_content("requirements.txt")

setup(
    name=geoplot.__name__,
    version=geoplot.__version__,
    description=geoplot.__description__,
    long_description=get_content("README.md"),
    long_description_content_type="text/markdown",
    author=geoplot.__author__,
    author_email=geoplot.__email__,
    maintainer=geoplot.__author__,
    maintainer_email=geoplot.__email__,
    classifiers=[el.lstrip() for el in """
        Development Status :: 3 - Alpha
        Intended Audience :: Science/Research
        Intended Audience :: System Administrators
        License :: OSI Approved :: MIT License
        Natural Language :: English
        Operating System :: OS Independent
        Programming Language :: Python
        Programming Language :: Python :: 3
        Programming Language :: Python :: 3.8
        Programming Language :: Python :: 3.9
        Topic :: Scientific/Engineering
    """.split('\n')],
    packages=find_packages(),
    zip_safe=False,
    install_requires=get_content("requirements.txt"),
    include_package_data=True)
