# tkinterGmshOGS-GINA_UpscaleMapping

## Dependencies

The code itself doesn't require any installation but it requires a Python installation in the computer system as well as the Python packages installed in the Python installation through 

```
pip install <package>
```

The list of package dependencies is as follows:

+ pathlib (pahtlib is now part of Python standard library as of Python v3.4) Thus it requires pip installation from v3.3 and earlier.
+ tkinter
+ time
+ scipy
+ pandas
+ numpy
+ ast
+ matplotlib
+ pykrige
+ sklearn
+ pyvista


## Description

The Code creates a Graphical User Interface (see GUI section below) for browsing the input files and selecting/typing the necessary input parameters. The spatial distribution of the input points across the mesh domain and under which elements they fall is illustrated in Figure 1.

![Figure 1 - Attribute upscaled and mesh .msh file nodes locations. The coloured points show the a locations and value of the property from the input .csv file whereas the black crosses show the mesh nodes. Every four nearby nodes form the edges of a 2-dimensional quadrilateral element.](https://github.com/benitez9rh/tkinterGmshOGS-GINA_UpscaleMapping/blob/main/GmshAttributeDist.png)

### Element value prediction

In the background, the code identifies the input points' X,Y locations and under which element the points will contribute to the prediction. This is dependent on the prediction method used.
The element value prediction depends on the averaging method chosen. If "Arithmetic Averaging" method is chosen, the points within the bounds of each particular element will be averaged and the result will be mapped to that element. If "Element Centre Kriging" is chosen, the input parameters used will form the major and minor continuity vectors which in turn make up the ellipse of correlation. The points within this ellipse will all contribute to the prediction of the point at the centre of the element with varying contributing weights; this may mean that points outside of the element or not all points inside the element will contribute to the prediction. The contributing weights of the points used and which points are used in this prediction both depend on the input parameters used.

![Figure 2 - Averaging method. In "Arithmetic Averaging", the points within the bounds of each particular element will be avergaed and the result will be mapped to that element. In the case of "Element Centre Kriging", the input parameters used will form the major and minor continuity vectors which in turn make up the ellipse of correlation. The points within this ellipse will all contribute to the prediction of the point at the centre of the element; this may mean that points outside of the element or not all points inside the element will contribute to the prediction, depending on the input parameters used.](https://github.com/benitez9rh/tkinterGmshOGS-GINA_UpscaleMapping/blob/main/PredictionMethod.png)

### Result

Once the predictions for each element and the mapping have been made, the spatial distribution of the property in the mesh will look like in Figure 3.

![Figure 3 - The mesh is explicit in this figure where the locations where the mesh lines cross correspond to the nodes. The colour correspond to the upscaled and mapped property from the input .csv file.](https://github.com/benitez9rh/tkinterGmshOGS-GINA_UpscaleMapping/blob/main/ModelAperture.png)


## Data

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


## GUI

Below, a snippet of the GUI is provided below.

![Figure 3 - The GUI for mapping a property or attribute onto a mesh.](https://github.com/benitez9rh/tkinterGmshOGS-GINA_UpscaleMapping/blob/main/GUI.PNG)

### Browse frame
The GUI allows you to browse for the .csv file containing the property data you wish to map onto the mesh as well as the mesh .msh file. Once you have done that, previews of the files are shwon in the corresponding preview entry boxes which allows for an easier input of the reading parameters without the need to open the files in third-party software. Check the columns you want to use from the file (containing the x and y locations as well as the property value) and type them in the "Columns to read" box separated by commas. If there is a header, type the header row value in the "Header row" entry box. IF there are any rows that you wish to skip, please also type the value in the "Rows to skip" entry box separated by commas if necessary. Leave blank if no rows to skip are necessary. Remember that Python indexing starts at 0 hence if you wish the first row to be read as the header, type 0 instead. 

### Prediction frame

In the "Prediction frame", you can chose the Material Group
The "Prediction method" refers to the algorithm used to calculate the element property value. As of version 1, only two are possible: "Arithmetic Averaging" and "Element Centre Kriging". If  "Arithmetic Averaging" is chosen, the "Variogram model frame" is hiden. If "Element Centre Kriging" is chosen, this option uses an PyKrige's Ordinary Kriging algorithm and takes the parameters from the "Variogram model frame" as input parameters for the algorithm.
The "Unsampled Element Attribute Value" refers to the value an element will take if there are no samples within its bounds. Options are "max", "min", "na" or a specific value.

### Variogram model frame

Ignore this frame if the "Element Averaging" prediction method is selected. The GUI should update automatically and disable these option anyway.
If a spatial continuity of the data has been performed a priori, you can use the spatial continuity parameters in the Ordinary Kriging algorithm prediction.
Otherwise, choose the "Variogram model" using the radio buttons. As of version 1, only the spherical model is available.
Choose the "major" and "minor continuity direction ranges", the "major continuity direction azimuth angle" measured clockwise from north, the "sill" and "nugget" parameters and type all these parameters in the respective entry boxes.


## Output

The code outputs a tab delimited .txt file, which can also be read as a .csv, containing two columns: the first column corresponds to the element number in the mesh and the second column corresponds to the mapped property value of each element. The file also contains a head which is readable by OGS. An example is shown below.

```
#MEDIUM_PROPERTIES_DISTRIBUTED
 $MSH_TYPE
  $LIQUID_FLOW
 $MMP_TYPE
  GEOMETRY_AREA
 $DIS_TYPE
  ELEMENT
 $DATA
  0  0.1
  1  0.2
  2  0.3
#STOP
// Averaging method: {prediction_method.get()}
```

The code can output the matplotlib.pyplot plots shown above.

## Pitfalls

The code testing performed thus far shows the code behaves well and runs quickly, especially when using "Arithmetic Averaging".
The code broke for lack of available memory during the PyKrige section when using "Element Centre Kriging" with more than 100,000 input points. This is due to the amount of memory necessary to build the kriging matrix inside the PyKrige.OrdinaryKriging code. Caution is advised when using large datasets.
