# Upscale and map scattered points onto OpenGeoSys GINA and Gmsh meshes formats: a Tkinter Graphical User Interface Python code

Gonçalo Benitez Cunha<sup>1</sup>

<sup>1</sup>The University of Edinburgh, Scotland, United Kingdom

Corresponding author: Gonçalo Cunha ([g.cunha@ed.ac.uk](g.cunha@ed.ac.uk))


## 1. Summary

OpenGeoSys (OGS) is an open-source numerical simulator for coupled thermal-hydraulic-mechanical-chemical (THMC) processes, particularly in porous and fractured media. OGS is used in several application modelling areas such as contaminant transport, regional and coastal hydrology, geothermal systems and energy storage, CO<sub>2</sub> sequestration and hydrogen storage and nuclear waste management and disposal (Kolditz et al. (2012), Bilke et al. (2022)). OGS utilises mesh files to discretise the domain for the numerical simulations, which are created by other software such as Gmsh (Geuzaine & Remacle (2020)), as well as element properties which are frequently heterogeneous across the domain.
The Graphical User Interface (GUI) presented in this paper offers the option to map a property onto a 2-dimensional quadrilateral mesh which OGS can then use to compute the numerical simulation, whilst offering some upscaling options including arithmetic averaging and spherical Ordinary Kriging through PyKrige (Murphy, (2014)). 

![tkinterGmshOGS-GINA_UpscaleMapping GUI](https://github.com/benitez9rh/tkinterGmshOGS-GINA_UpscaleMapping/blob/main/GUI.PNG)

*Figure 1 - tkinter GUI for property upscale and mapping onto Gmsh and  OGS GINA's mesh formats*

## 2. Statement of need

Finite Element Method (FEM) numerical models rely on media properties, including in coupled thermal-hydraulic-mechanical-chemical (THMC) processes. Frequently, these properties are not homogeneous across the domain. In these situations, there is often a need to map a property from a .csv point cloud format onto the model's mesh for numerical simulation. This would look like a property value per element .txt or .csv file, requiring prior data structure formats and coding knowledge.
In this paper, we try to partially bridge that gap by presenting a user-friendly graphical user interface (GUI) writen in Python using tkinter. The code maps a property onto a 2-dimensional quadrilateral element mesh whilst providing some upscaling options and can handle OpenGeoSys GINA's and Gmsh formats.
The GUI streamlines the experience to browse the property and mesh files making it simple and easy to use and quickly providing the user with a file to input in OGS to continue with the numerical simulation without the need to code.

