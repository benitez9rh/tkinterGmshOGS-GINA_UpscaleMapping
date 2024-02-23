---
title: 'Upscale and map scattered points onto OpenGeoSys GINA and Gmsh meshes formats: a Tkinter Graphical User Interface Python code'

tags:
- Python
- Tkinter GUI
- OpenGeoSys
- Gmsh
- Finite Element Method
- Numerical modelling

authors:
- name: Gonçalo Benitez Cunha
- corresponding: true
- orcid: 0009-0007-5441-3790
- affiliation: 1

affiliations:
- name: The University of Edinburgh, Scotland, United Kingdom
  index: 1
  
date: 15 February 2024

bibliography: paper.bib

---


## Summary

OpenGeoSys (OGS) [@Kolditz2012] is an open-source numerical simulator for coupled thermal-hydraulic-mechanical-chemical (THMC) processes, particularly in porous and fractured media. OGS is used in several application modelling areas such as contaminant transport, regional and coastal hydrology, geothermal systems and energy storage, CO<sub>2</sub> sequestration and hydrogen storage and nuclear waste management and disposal [@Kolditz2012; @Bilke2022]. OGS utilises mesh files to discretise the domain for the numerical simulations, which are created by other software such as Gmsh [@Gmsh2020], as well as element properties which are frequently heterogeneous across the domain.
This paper presents a Graphicalk User Interface (GUI), coded in tkinter [@Lundh1999], offers the option to map a property onto a 2-dimensional quadrilateral mesh which OGS can then use to compute the numerical simulation, whilst offering some upscaling options. These options include arithmetic averaging and Ordinary Kriging using a spherical semi-variogram model through PyKrige [@Murphy2022].


## Statement of need

Numerical modelling techniques, such as Finite Element Method (FEM), rely on media properties, including in coupled thermal-hydraulic-mechanical-chemical (THMC) processes. Frequently, these properties are not homogeneous across the domain. In these situations, there is often a need to map a property from a .csv point cloud format onto the model's mesh for numerical simulation. In other words, a .txt or .csv file containing the upscaled property value per element, which would require prior data structure formats and coding knowledge to construct.
In this paper, we try to partially bridge that gap by presenting a user-friendly graphical user interface (GUI) writen in Python using tkinter [@Lundh1999]. The code maps a property onto a 2-dimensional quadrilateral element mesh whilst providing some upscaling options and can handle OpenGeoSys's GINA [@Kolditz2016] and Gmsh [@Gmsh2020] .msh formats.
The GUI streamlines the experience to browse the property and mesh files making it simple and easy to use and quickly providing the user with a file to input in OGS to continue with the numerical simulation without the need to code.


## Acknowledgements

This work was possible due to the sponsorship of Quintessa, the School of Geosciences at the The University of Edinburgh and Nuclear Waste Services, part of the Nuclear Decommissioning Authority in the United Kingdom.
The code uses various other Python's [@VanRossum1995] packages such as PyKrige [@Murphy2022], Scipy [@Virtanen2020], Pandas [@McKinney2010], Numpy [@Harris2020] and matplotlib [@Hunter2007]. Other smaller and non-published code contributions used are cited and acknowledged accordingly as comments and/or links in the code.


## References


