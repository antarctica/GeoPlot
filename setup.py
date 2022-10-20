from importlib.metadata import entry_points
from setuptools import setup, find_packages
import glob

import bas_geoplot


def get_content(filename):
    with open(filename, "r") as fh:
        return fh.read()


requirements = get_content("requirements.txt")

setup(
    name=bas_geoplot.__name__,
    version=bas_geoplot.__version__,
    description=bas_geoplot.__description__,
    long_description=get_content("README.md"),
    long_description_content_type="text/markdown",
    author=bas_geoplot.__author__,
    author_email=bas_geoplot.__email__,
    maintainer=bas_geoplot.__author__,
    maintainer_email=bas_geoplot.__email__,
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
    data_files= glob.glob('bas_geoplot/config/**'),
    entry_points={
        'console_scripts': [
            "plot_mesh=bas_geoplot.cli:plot_mesh_cli"
        ]
    },
    packages=find_packages() + ["bas_geoplot.config"],
    zip_safe=False,
    install_requires=requirements,
    include_package_data=True)
