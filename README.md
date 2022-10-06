![](logo.jpg)
# GeoPlot
>GeoPlot is an interactive and static plotting toolkit that is leveraged by several projects across the British Antarctic Survey Artificial Intelligence Lab (BAS AI) used in combinations with additional software packages being developed. 

## Installation
The software package can be installed by downloading the github repo and running
```
python setup.py install
```
or equally installed directly using
```
pip install geoplot
```

## Examples
There are two types of plotting functions included in the software package: `Static` and `Interactive`. Outlined below are examples how each of these sections of this software package can be used, with further documentation given in the doc strings of the classes. For both classes the plotting is run by applying consecutive layers to plot by consecutive function calls. For each of the software packages being generated within the BAS AI Lab many of the datasets have a list of plotting standard, be it colour for a plotted mesh or size of a marker, so we leverage a set of standard as defined within the configuration files for each plotting type, given in the `./geoplot/config` folder.
### *Interactive*
The interactive plotting functions leverage the python `folium` package to generate interactive `html` files that can be generated locally and run using a web browser. 
``` python
from geoplot.interactive import iMap
mp = iMap(title='Example Plot')
mp.Maps(mesh,'SIC',predefined='SIC')
mp.Points(itinary,'Waypoints')
mp.MeshInfo(mesh,'Mesh Info',show=False)
mp.show()
```
In this example the `mp = iMap(title='Example Plot')` generates the initial basemap for a region of interest, with subseqent lines of code adding layers on-to this basemap. `predefined` is used to load predefined in the configuration file.

### *Static*
The static plotting functions leverage the python `cartopy` package to generate static files. For simplicity we follow simialar structure to the interactive plots, with an example coding section given by
``` python
from geoplot.static import sMap
mp = sMap(title='Example Plot')
mp.Maps(mesh,predefined='SIC')
mp.Points(itinary)
mp.savefig('Figure.pdf')
```

## Development & Contributions
Development of software package is conducted by the [BAS AI Lab](https://www.bas.ac.uk/team/science-teams/ai-lab/).

For contruibutions and feature additions please contact [jonsmi@bas.ac.uk](jonsmi@bas.ac.uk).

Distributed under the MIT license. See ``LICENSE`` for more information.


[version]: https://img.shields.io/GeoPlot/v/datadog-metrics.svg?style=flat-square
[downloads]: https://img.shields.io/GeoPlot/dm/datadog-metrics.svg?style=flat-square
