# Upscale and map scattered points onto OpenGeoSys GINA and Gmsh meshes formats: a Tkinter Graphical User Interface Python code

Gonçalo Benitez Cunha<sup>1</sup>

<sup>1</sup>The University of Edinburgh, Scotland, United Kingdom

Corresponding author: Gonçalo Cunha ([g.cunha@ed.ac.uk](g.cunha@ed.ac.uk))


## 1. Summary

OpenGeoSys (OGS) is an open-source numerical simulator for coupled thermal-hydraulic-mechanical-chemical (THMC) processes, particularly in porous and fractured media. OGS is used in several application areas such as contaminant transport, regional and coastal hydrology, geothermal systems and energy storage, CO<sub>2</sub> sequestration and hydrogen storage and nuclear waste management and disposal (Kolditz et al. (2012), Bilke et al. (2022)). OGS utilises mesh files to discretise the domain for the numerical simulations, which are created by other software such as Gmsh (Geuzaine & Remacle (2020)), as well as element properties which are frequently heterogeneous across the domain. The Graphical User Interface (GUI) presented in this paper offers the option to map  a property onto a 2-dimensional quadrilateral mesh which OGS can then use to compute the numerical simulation, whilst offering some upscaling options including arithmetic averaging and spherical Ordinary Kriging through pykrige (Murphy, (2014)). 

## 2. Statement of need

Finite Element Method (FEM) numerical models rely on media properties, including in coupled thermal-hydraulic-mechanical-chemical (THMC) processes. Frequently, these properties are not homogeneous across the domain. In these situations, there is often a need to map a property from a .csv point cloud format onto the model's mesh for numerical simulation. This would look like a property value per element .txt or .csv file, requiring prior data structure formats and coding knowledge.
In this paper, we try to partially bridge that gap by presenting a user-friendly graphical user interface (GUI) writen in Python using tkinter. The code maps a property onto a 2-dimensional quadrilateral element mesh whilst providing some upscaling options and handles OpenGeoSys GINA's and Gmsh formats. The GUI streamlines the experience to browse the property and mesh files making it simple and easy to use.

## 3. Graphical User Interface



## 4. Data

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


