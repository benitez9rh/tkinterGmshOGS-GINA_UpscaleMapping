# Upscale and map scattered points onto OpenGeoSys GINA and Gmsh meshes formats: a Tkinter Graphical User Interface Python code

Gonçalo Benitez Cunha<sup>1</sup>

<sup>1</sup>The University of Edinburgh, Scotland, United Kingdom

Corresponding author: Gonçalo Cunha ([g.cunha@ed.ac.uk](g.cunha@ed.ac.uk))


## 1. Summary

OpenGeoSys (OGS) is an open source numerical simulator for coupled thermo-hydraulic-mechanical-chemical (THMC) processes, particularly in porous and fractured media. Due to the fact that it is open source, it requires prior data structure formats and coding knowledge, apart from the scientific knowledge intrinsic to THMC coupled processes and numerical modelling. 
In this paper, we present a graphical user interface writen in Python using tkinter

## 2. Statement of need

OpenGeoSys (OGS) is an open source numerical simulator for coupled thermo-hydraulic-mechanical-chemical (THMC) processes in porous and fractured media. OGS is used in several application areas such as contaminant transport, regional and coastal hydrology, geothermal systems and energy storage, CO<sub>2</sub> sequestration and hydrogen storage and nuclear waste management and disposal (Kolditz et al. (2012), Bilke et al. (2022)).


## 3. Data

The data required needs to be in a format that is readable by thePython Pandas package, such as tab delimited .txt or .csv files, and they must contain the required XYZ spatial information. The input mesh (.msh) files from Gmsh and imported (re-formatted) by GINA are also shown below.
Below an example of the format required:

```
***  .csv  ***
x	y	aperture	avg	var	std
-22.5	-32.5	0.49977896	0.515901725	0.005059804	0.071132302
-22.5	-31.5	0.48803121	0.508111098	0.004592228	0.067765978         
...     ...     ...         ...         ...         ...

***  GINA's .msh  ***
$NODES
  8186 
0 -22.5 -32.5 1
1 -22.5 -31.5 1
...
$ENDNODES
$ELEMENTS
  8184 
0 0 tri 0 3977 4220
1 0 tri 4220 3977 4465
...
$ENDELEMENTS

***  Gmsh's .msh  ***
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
10440
1 -40.1906200200441 107.189905102499 0
...
10440 70.17557769416671 106.1857036355631 0
$EndNodes
$Elements
10829
1 15 2 0 0 1
...
10829 3 2 0 3 10440 493 8 494
$EndElements
```


