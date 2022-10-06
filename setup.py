from setuptools import setup, find_packages

import geo_plot


def get_content(filename):
    with open(filename, "r") as fh:
        return fh.read()


requirements = get_content("requirements.txt")

setup(
    name=geo_plot.__name__,
    version=geo_plot.__version__,
    description=geo_plot.__description__,
    long_description=get_content("README.md"),
    long_description_content_type="text/markdown",
    author=geo_plot.__author__,
    author_email=geo_plot.__email__,
    maintainer=geo_plot.__author__,
    maintainer_email=geo_plot.__email__,
    url="https://github.com/antarctica/GeoPlot",
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
    packages=find_packages() + ["geo_plot.config"],
    zip_safe=False,
    install_requires=get_content("requirements.txt"),
    include_package_data=True)
