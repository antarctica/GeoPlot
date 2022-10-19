![](logo.jpg)
<!-- ![Version][version] ![Downloads][downloads] -->

>GeoPlot is an interactive and static plotting toolkit that is leveraged by several projects across the British Antarctic Survey Artificial Intelligence Lab (BAS AI) used in combinations with additional software packages being developed. 

## Manual Pages

All information on the installation and usage of this software can be found at the manual pages [LINK](https://github.com/antarctica/GeoPlot/wiki).

## Installation
The software package can be installed by downloading the github repo by running the examples below. This software package requires GDAL and Fiona to be pre-installed, this is not an issue with MacOS/Linux installs but additional steps must be taken for Windows install. These additional steps are as follows

windows only:
```
    pip install pipwin
    pipwin install gdal
    pipwin install fiona
    pipwin install cartopy
```

For the installation of the GeoPlot run one of the following:

source install:
```
python setup.py install
```

pip install:
```
pip install bas-geoplot
```


## Development & Contributions
Development of software package is conducted by the [BAS AI Lab](https://www.bas.ac.uk/team/science-teams/ai-lab/).

For contruibutions and feature additions please contact [jonsmi@bas.ac.uk](jonsmi@bas.ac.uk).


## License
Distributed under the MIT license. See ``LICENSE`` for more information.

[version]: https://img.shields.io/GeoPlot/v/datadog-metrics.svg?style=flat-square
[downloads]: https://img.shields.io/GeoPlot/dm/datadog-metrics.svg?style=flat-square