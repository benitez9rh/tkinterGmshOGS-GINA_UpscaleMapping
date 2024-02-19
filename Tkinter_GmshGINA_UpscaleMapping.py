# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 16:37:50 2024

@author: s2132627
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:28:20 2023
Created on Mon May  8 10:32:54 2023
Created on Wed Apr 26 11:46:02 2023
@author: s2132627
"""


"""
This script maps initial point (x,y,z, aperture) data to mesh data (nodes, elements) so that we can create apertures and hydraulic conductivities (currently hardcoded for the cubic law) text files for the mesh elements (as opposed to spatial points) that OGS can read.

In this version, I use the elements' central points in the .msh file and krig the aperture from the .csv file. It requires therefore the spatial continuity parameters for the kriging algorithm. Alternatively, as with the last version, we can use the element's bounds x,y,z to filter all the .csv points that fall within that element and average them (works only for quad elements at the moment because I have been having issues with the winding number's algorithm). This way I don't have to faff about with KDtrees and such. Much simpler, probably more accurate.




Below are an example of the .csv file and the (GINA's MeshFormat 2.2) .msh and a Gmsh .msh file sfor reference.
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
 $ELEMENTS
  8184 
0 0 tri 0 3977 4220
1 0 tri 4220 3977 4465

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
697 3 2 0 1 597 598 686 685
...
10829 3 2 0 3 10440 493 8 494
$EndElements
"""
import os
import pathlib
import sys
from tkinter import *               #Import everything from Tkinter
from tkinter import ttk
import tkinter.messagebox           #Import messagebox method
from tkinter import filedialog
from time import time, ctime
import time
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import pandas as pd
from numpy import *
import numpy as np
from ast import literal_eval
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from pykrige.ok import OrdinaryKriging
import pykrige.kriging_tools as kt
# import threading
# import subprocess
# from scipy import stats
# from scipy.stats import zscore # imports the normal score method used to get rid of the outliers in the data
# from scipy.stats import lognorm, norm
# import csv
# import math as m
# from sklearn.metrics import r2_score
# import pyvista as pv
# from pyvista import examples
# from Euler3DRotation import rotate_np


'''User Defined Global Variables'''
output_extension = ".txt"
eamthreshold = 1E-30
save = True

"""System Variables"""
cmap = plt.cm.plasma                    # color map
extension_xyz = ".xyz"
extension_csv = ".csv"
extension_png = ".png"
extension_txt = ".txt"

'''Tkinter Variables'''
root = Tk()
BROWSE_frame = LabelFrame(root, text = "Browse Frame")
BROWSE_frame.grid(row = 1, column = 1, columnspan = 3, sticky = EW, pady = 10, padx = 10)
CSV_frame1 = LabelFrame(root, highlightbackground="black", highlightthickness=1)
CSV_frame1.grid(row = 3, column = 1, pady = 20, padx = 10)
MSH_frame = Frame(root, highlightbackground="black", highlightthickness=1)
MSH_frame.grid(row = 3, column = 3, pady = 20, padx = 10)
PRED_frame = LabelFrame(root, text = "Prediction Frame")
PRED_frame.grid(row = 2, column = 1, pady = 20, padx = 10)
VARMODEL_frame = LabelFrame(root, text = "Variogram Model Frame (ignore if Element Averaging is selected in the Prediction Frame)")
VARMODEL_frame.grid(row = 2, column = 2, pady = 20, padx = 10)
prediction_method = StringVar(value = "Element Centre Kriging")
conductivity_method = StringVar(value = "Cubic Law")
varmdl_method = StringVar(value = "Spherical")

'''Methods'''
def simplest_type(s): # ast 's literal_eval converts numericals to their appropriate type. Say you have "2" and "3.14" strings, they will be converted to int and float respectively. The issue is that it breaks when it sees an actual string, thus this simplest_type function. # From https://stackoverflow.com/questions/60213604/automatically-convert-string-to-appropriate-type
    try:
        return literal_eval(s)
    except:
        return s
def normal_round(n):
    if n - m.floor(n) < 0.5:
        return m.floor(n)
    return m.ceil(n)
def prec(df):
    precision = lambda x: (len(x)-1, len(x.split('.')[1]))
    precision = df.astype(str).applymap(precision).max().max()[1]
    #print(precision)
    return precision
def nsmall(a, n):
    return np.partition(a, n)[n]                                                #finds the n-th smallest value
def msh_of_SplitFunction_v4(path, input_file_name):
    """
    Initially built to create a pandas Data Frame from the .msh_of file from GINA's OGS pre-processor (hence the name), it now is able to distinguish between that and a Gmsh .msh file to create the same pandas df.
    
    Parameters
    ----------
    path : String
        system path address to the file.
    input_file_name : TYPE
        system file name provided without the extension. In the case of GINA's .msh_of, please provide the .msh file and the code will look for the .msh_of file itself.

    Returns
    -------
    dfN : Pandas Data Frame
        Pandas Data Frame with the nodes and their xy(z) locations.
    dfE : Pandas Data Frame
        Pandas Data Frame with the elements and their corresponding nodes.
    ENC_df : Pandas Data Frame
        Pandas Data Frame with the elements and their values.
    MGs : Pandas Data Frame
        Pandas Data Frame with the material groups of each element.
    Etype : String
        The type of element.

    """
    # #############################  Get the Nodes and Elements from the .msh (GINA's) file without having to export and read back .txt files #############################
    
    with open(path + input_file_name + ".msh", 'r') as file:
        
        # Check the file format
        head = file.readlines()[0:50] # Check the first 50 lines
        print(head)
        Format = "GINA"
        for l in head:    
            if "$MeshFormat" in l:
                Format = "Gmsh"
        if Format == "GINA":
            words = ["$NODES", "$ELEMENTS"]
        elif Format == "Gmsh":
            words = ["$Nodes", "$Elements"]
        print(Format)
        
    with open(path + input_file_name + ".msh", 'r') as file:                            # The below code block doesn't run it I put it in the same with open block above. I don't know why
        # This Loop is not working in tkinter for whatever reason so I hardcoded it
        # for word in words: 
        #     exec(f"{word[word.rfind('$')+1:]} = False") # Creates a False flag for each word in words
        NODES = False
        ELEMENTS = False
        for line in file:
            # print(line)
            casematch = next((word for word in words if word in line), False)
            if casematch:
                # for word in words: # Creates a False flag for each word in words  
                #     exec(f"{word[word.rfind('$')+1:]} = False")
                NODES = False
                ELEMENTS = False
                # exec(f"{match[match.rfind('$')+1:]} = True") # Sets the match word flag to True
                print(casematch)
                print(NODES)
                print(ELEMENTS)
                # if casematch.find("NODES") != -1 or casematch.find("Nodes") != -1: # If line contains "NODES"
                if "nodes" in casematch.lower():
                    NODES = True
                    Nnodes = int(next(file))
                    print(Nnodes)
                    Ns = []
                # elif casematch.find("ELEMENTS") != -1 or casematch.find("Elements") != -1: # If line contains "ELEMENTS"
                elif "elements" in casematch.lower():
                    ELEMENTS = True
                    Nelements = int(next(file))
                    print(Nelements)
                    Es = []
            else:
                if Format == "GINA":
                    if line.find("#STOP") != -1:
                        break
                    elif NODES == True:
                        Ns.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                    elif ELEMENTS == True:
                        Es.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                elif Format == "Gmsh":
                    if "$End" in line: # if line contains
                        NODES = False
                        ELEMENTS = False
                    if NODES == True:
                        Ns.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
                    elif ELEMENTS == True:
                        Es.append(  [simplest_type(v) for v in line[:line.rfind("\n")].split(" ") ] )
        
    dfN = pd.DataFrame(Ns, columns = ["NodeTag","x","y","z"]) # reads only the nodes
    # dfN = dfN.round(prec(df)) # This is necessary because of precision changes, later on when I am trying to find the aperture corresponding to the x,y,z values of a node, those x,y,z won't match.
    if Format == "GINA":
        Etype = Es[0][2] # Get E type (tri, quad)
        EsMaxLen = max([len(i) for i in Es])
        names = ["ElementNumber", "MaterialGroup", "ElementType"] 
        [names.append (f"Node{i+1}") for i in range(EsMaxLen - len(names))] # Add entries to the names list corresponding to the unknown number of nodes
        dfE = pd.DataFrame(Es, columns = names); #dfE = dfE.round(prec(df))
        EmaxN = max([(simplest_type(n[n.rfind("Node")-1:])) for n in list(dfE.columns) if n.find("Node") != -1]) # Count the number of columns with "Node"
        print(dfN)
        print(dfE)
        if os.path.isfile(path + input_file_name + ".msh_of"):
            with open(path + input_file_name + ".msh_of") as f:
                contents = f.read()
            f.close()
            Esplit = contents.split("\n\n")
            Esplit = [ s.split("\n") for s in Esplit] # Split contents into elements
            ENC = [[simplest_type(l2.strip()) for l2 in l1] for l1 in Esplit] # Elements' strings cleaned and converted into numbers    # Instead of for i, l1 in enumerate(f2split):
            EsMaxLen = max([len(i) for i in ENC])
            EsC_columns = ["MaterialGroup", "ElementNumber"] 
            [EsC_columns.append (f"Node{i+1}") for i in range(EsMaxLen - len(EsC_columns) - 5)] # Add entries to the names list corresponding to the unknown number of nodes. To get the unknown number of nodes, loop through a range to N where N is the number of total values read from the file minus the number of values that are known to make up the file.
            EsC_columns.extend (["ElementCentreX-Coordinate", "ElementCentreY-Coordinate", "ElementCentreZ-Coordinate", "Dunno", "Dunno2" ]) # Elements' Centres
            if "qua" in Etype.lower():
                ENC_df = pd.DataFrame(ENC[1:], columns = EsC_columns) # Start at row1 ignoring row 0 because the latter is not an Element
                # ENC_df_ElementCentrexCoordinate_round = round(ENC_df["Element Centre x-Coordinate"])
                ENC_df["GridX-coordinate"]  = (round(ENC_df["ElementCentreX-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))                                        # Round because there might be a precision issue. This rounding does not change the initial dataframe. -1 because the ranking function starts at 1
                ENC_df["GridX-coordinate"]  = ENC_df["GridX-coordinate"] / (np.sort(ENC_df["GridX-coordinate"].unique())[1] - np.sort(ENC_df["GridX-coordinate"].unique())[0])      # rank function assigns the same rank to duplicates but jumps that amount of duplicates when the next value is different, which we don't want. Thus we are calculating the difference between sequential values and dividing..
                ENC_df["GridY-coordinate"]  = (round(ENC_df["ElementCentreY-Coordinate"]).rank(method = "min") - 1).apply(lambda x: floor(x))
                ENC_df["GridY-coordinate"]  = ENC_df["GridY-coordinate"] / (np.sort(ENC_df["GridY-coordinate"].unique())[1] - np.sort(ENC_df["GridY-coordinate"].unique())[0])
        else:      
            ENC_df = dfE
            ENC_df["ElementCentreX-Coordinate"] = 0
            ENC_df["ElementCentreY-Coordinate"] = 0
            ENC_df["ElementCentreZ-Coordinate"] = 0
            ENC_df["GridX-coordinate"] = 0
            ENC_df["GridY-coordinate"] = 0
            for row in ENC_df.itertuples():
# =============================================================================
#                 # Deprecated because below I don't need to hardcode
#                 if "tri" in Etype.lower():
#                     Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3") ]
#                 if "qua" in Etype.lower():
#                     Ns = [ getattr(row, "Node1"), getattr(row, "Node2"), getattr(row, "Node3"), getattr(row, "Node4") ]
# =============================================================================
                Ns = [getattr(row, f"Node{i+1}") for i in range(EmaxN)]
                Ns = dfN[dfN["NodeTag"].isin(Ns)]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreX-Coordinate']  = Ns.x.describe()["mean"]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreY-Coordinate']  = Ns.y.describe()["mean"]
                ENC_df.at[getattr(row, "Index"), 'ElementCentreZ-Coordinate']  = Ns.z.describe()["mean"]
            ENC_df["GridX-coordinate"]  = np.argsort(ENC_df["ElementCentreX-Coordinate"])
            ENC_df["GridY-coordinate"]  = np.argsort(ENC_df["ElementCentreY-Coordinate"])
            ENC_df["GridZ-coordinate"]  = np.argsort(ENC_df["ElementCentreZ-Coordinate"])
        MGs = np.sort(ENC_df["MaterialGroup"].unique() )# np.array containing all used Material Groups' numbers
    elif Format == "Gmsh":
        EsMaxLen = max([len(i) for i in Es])
        names = ["ElementNumber", "ElementTypeNumber", "ElementUnknown1", "ElementUnknown2", "MaterialGroup"] # Element number, Element type, Element unknown characteristic 1, Element unknown characteristic 2, Element Material Group
        [names.append (f"Node{i+1}") for i in range(EsMaxLen - len(names))] # Add entries to the names list corresponding to the unknown number of nodes
        dfE = pd.DataFrame(Es, columns = names)
        
        Gmshtypes_unique = dfE["ElementTypeNumber"].unique().tolist() # The second column has the element identifier. See https://gitlab.onelab.info/gmsh/gmsh/blob/master/src/common/GmshDefines.h as of 12.02.2024
        Gmshtypes = [
                    ['MSH_LIN_2', 1],         ['MSH_TRI_3', 2],         ['MSH_QUA_4', 3],         ['MSH_TET_4', 4],         ['MSH_HEX_8', 5],         ['MSH_PRI_6', 6],         ['MSH_PYR_5', 7],         ['MSH_LIN_3', 8],
                    ['MSH_TRI_6', 9],         ['MSH_QUA_9', 10],         ['MSH_TET_10', 11],         ['MSH_HEX_27', 12],         ['MSH_PRI_18', 13],         ['MSH_PYR_14', 14],         ['MSH_PNT', 15],         ['MSH_QUA_8', 16],
                    ['MSH_HEX_20', 17],         ['MSH_PRI_15', 18],         ['MSH_PYR_13', 19],         ['MSH_TRI_9', 20],         ['MSH_TRI_10', 21],         ['MSH_TRI_12', 22],         ['MSH_TRI_15', 23],         ['MSH_TRI_15I', 24],
                    ['MSH_TRI_21', 25],         ['MSH_LIN_4', 26],         ['MSH_LIN_5', 27],         ['MSH_LIN_6', 28],         ['MSH_TET_20', 29],         ['MSH_TET_35', 30],         ['MSH_TET_56', 31],         ['MSH_TET_22', 32],
                    ['MSH_TET_28', 33],         ['MSH_POLYG_', 34],         ['MSH_POLYH_', 35],         ['MSH_QUA_16', 36],         ['MSH_QUA_25', 37],         ['MSH_QUA_36', 38],         ['MSH_QUA_12', 39],         ['MSH_QUA_16I', 40],
                    ['MSH_QUA_20', 41],         ['MSH_TRI_28', 42],         ['MSH_TRI_36', 43],         ['MSH_TRI_45', 44],         ['MSH_TRI_55', 45],         ['MSH_TRI_66', 46],         ['MSH_QUA_49', 47],         ['MSH_QUA_64', 48],
                    ['MSH_QUA_81', 49],         ['MSH_QUA_100', 50],         ['MSH_QUA_121', 51],         ['MSH_TRI_18', 52],         ['MSH_TRI_21I', 53],         ['MSH_TRI_24', 54],         ['MSH_TRI_27', 55],         ['MSH_TRI_30', 56],
                    ['MSH_QUA_24', 57],         ['MSH_QUA_28', 58],         ['MSH_QUA_32', 59],         ['MSH_QUA_36I', 60],         ['MSH_QUA_40', 61],         ['MSH_LIN_7', 62],         ['MSH_LIN_8', 63],         ['MSH_LIN_9', 64],
                    ['MSH_LIN_10', 65],         ['MSH_LIN_11', 66],         ['MSH_LIN_B', 67],         ['MSH_TRI_B', 68],         ['MSH_POLYG_B', 69],         ['MSH_LIN_C', 70],         ['MSH_TET_84', 71],         ['MSH_TET_120', 72],
                    ['MSH_TET_165', 73],         ['MSH_TET_220', 74],         ['MSH_TET_286', 75],         ['MSH_TET_34', 79],         ['MSH_TET_40', 80],         ['MSH_TET_46', 81],         ['MSH_TET_52', 82],         ['MSH_TET_58', 83],
                    ['MSH_LIN_1', 84],         ['MSH_TRI_1', 85],         ['MSH_QUA_1', 86],         ['MSH_TET_1', 87],         ['MSH_HEX_1', 88],         ['MSH_PRI_1', 89],         ['MSH_PRI_40', 90],         ['MSH_PRI_75', 91],
                    ['MSH_HEX_64', 92],         ['MSH_HEX_125', 93],         ['MSH_HEX_216', 94],         ['MSH_HEX_343', 95],         ['MSH_HEX_512', 96],         ['MSH_HEX_729', 97],         ['MSH_HEX_1000', 98],         ['MSH_HEX_32', 99],
                    ['MSH_HEX_44', 100],         ['MSH_HEX_56', 101],         ['MSH_HEX_68', 102],         ['MSH_HEX_80', 103],         ['MSH_HEX_92', 104],         ['MSH_HEX_104', 105],         ['MSH_PRI_126', 106],         ['MSH_PRI_196', 107],
                    ['MSH_PRI_288', 108],         ['MSH_PRI_405', 109],         ['MSH_PRI_550', 110],         ['MSH_PRI_24', 111],         ['MSH_PRI_33', 112],         ['MSH_PRI_42', 113],         ['MSH_PRI_51', 114],         ['MSH_PRI_60', 115],
                    ['MSH_PRI_69', 116],         ['MSH_PRI_78', 117],         ['MSH_PYR_30', 118],         ['MSH_PYR_55', 119],         ['MSH_PYR_91', 120],         ['MSH_PYR_140', 121],         ['MSH_PYR_204', 122],         ['MSH_PYR_285', 123],
                    ['MSH_PYR_385', 124],         ['MSH_PYR_21', 125],         ['MSH_PYR_29', 126],         ['MSH_PYR_37', 127],         ['MSH_PYR_45', 128],         ['MSH_PYR_53', 129],         ['MSH_PYR_61', 130],         ['MSH_PYR_69', 131],
                    ['MSH_PYR_1', 132],         ['MSH_PNT_SUB', 133],         ['MSH_LIN_SUB', 134],         ['MSH_TRI_SUB', 135],         ['MSH_TET_SUB', 136],         ['MSH_TET_16', 137],         ['MSH_TRI_MINI', 138],         ['MSH_TET_MINI', 139],
                    ['MSH_TRIH_4', 140]
                    ] # Added the corresponding number of nodes in  some of the mesh element types.
        Gmshtypes_df = pd.DataFrame(Gmshtypes, columns = ["ElementType" , "ElementTypeNumber"])
        # Gmshtypes_df['No_Nodes'] = Gmshtypes_df['No_Nodes'].fillna(0).astype(int) # Not absolutely necessary. The objective is to convert theh number of nodes to int datatype but because most types don't have a number of nodes associated, pandas converts the whole column to float and those rows that don't have a node value become NaNs.
        for i in Gmshtypes_unique:
            exec ( f"Es{i} = dfE[dfE['ElementTypeNumber'] == {i}]" )  # Separate Gmsh elements (all points, lines, any time of element) into different pd dfs
        dfE = dfE[dfE["ElementTypeNumber"].isin([15, 1]) == False]    # Subset of the DataFrame excluding all mesh points (dfE["EtypeNo"] = 15) and mesh lines(dfE["EtypeNo"] = 1)
        ENC_df = dfE
        ENC_df["ElementCentreX-Coordinate"] = 0
        ENC_df["ElementCentreY-Coordinate"] = 0
        ENC_df["ElementCentreZ-Coordinate"] = 0
        ENC_df["GridX-coordinate"] = 0
        ENC_df["GridY-coordinate"] = 0
        for row in ENC_df.itertuples():
            Ns = [ getattr(row, f"Node{n+1}") for n in range( len([n for n in ENC_df.columns if n.find("Node") != -1]) )]       # Create list with the nodes' numbers of each element by looping through all the column names that contain "Node". This way we don't need to know a priori the amount of nodes and the type of element
            Ns = dfN[dfN["NodeTag"].isin(Ns)]                                                                                   # Creates a pd df with only the nodes from Ns list above
            ENC_df.at[getattr(row, "Index"), 'ElementCentreX-Coordinate']  = Ns.x.describe()["mean"]                          #   Calculates the mean to find the centre
            ENC_df.at[getattr(row, "Index"), 'ElementCentreY-Coordinate']  = Ns.y.describe()["mean"]
            ENC_df.at[getattr(row, "Index"), 'ElementCentreZ-Coordinate']  = Ns.z.describe()["mean"]
        ENC_df["GridX-coordinate"]  = np.argsort(ENC_df["ElementCentreX-Coordinate"])                                        # Sorts the elements by grid coordinates
        ENC_df["GridY-coordinate"]  = np.argsort(ENC_df["ElementCentreY-Coordinate"])
        MGs = np.sort(ENC_df["MaterialGroup"].unique() )# np.array containing all used Material Groups' numbers                # Creates a np array with all the Material Groups and shorts it
        Etype = Gmshtypes_df.loc[Gmshtypes_df['ElementTypeNumber'] == dfE["ElementTypeNumber"].unique()[0], 'ElementType'].iloc[0]
        
    return dfN, dfE, ENC_df, MGs, Etype
class MidpointNormalize(colors.Normalize):
	"""
	Normalise the colorbar so that diverging bars work there way either side from a prescribed midpoint value)

	e.g. im=ax1.imshow(array, norm=MidpointNormalize(midpoint=0.,vmin=-100, vmax=100))
    # set the colormap and centre the colorbar. Thanks to http://chris35wills.github.io/matplotlib_diverging_colorbar/
	"""
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)
	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
class Spinner:
    """ Thanks to Victor Moyseenko https://stackoverflow.com/a/39504463    """
    busy = False
    delay = 0.1
    @staticmethod
    def spinning_cursor():
        while 1: 
            for cursor in '|/─\\': yield cursor
    def __init__(self, delay=None):
        self.spinner_generator = self.spinning_cursor()
        if delay and float(delay): self.delay = delay
    def spinner_task(self):
        while self.busy:
            sys.stdout.write(next(self.spinner_generator))
            sys.stdout.flush()
            time.sleep(self.delay)
            sys.stdout.write('\b')
            sys.stdout.flush()
    def __enter__(self):
        self.busy = True
        threading.Thread(target=self.spinner_task).start()
    def __exit__(self, exception, value, tb):
        self.busy = False
        time.sleep(self.delay)
        if exception is not None:
            return False

'''Tkinter Methods'''
def prediction_method_selection():
    prediction_method_text = f"**{str(prediction_method.get())}** Prediction method selected"
    prediction_method_label.config(text = prediction_method_text)
    if prediction_method.get() == "Element Centre Kriging":
        VARMODEL_frame.grid(row = 2, column = 2, pady = 20, padx = 10)
    else:
        VARMODEL_frame.grid_forget()
def varmdl_method_selection():
    varmdl_method_text = f"**{str(varmdl_method.get())}** Variogram Model selected"
    varmdl_method_label.config(text = varmdl_method_text)
def description():
    tkinter.messagebox.showinfo(title = "Description", message = " Simply use the 'Browse' button to browse to the output .CSV file from Gmsh.\n\n Then click 'Format', or alternatively 'File' > 'Format' and a new file will be created in the same directory formatted and ready to input in GINA.")
def browseCSV1():
    root.CSVinputfilepath1 = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".csv files","*.csv"), (".txt files", "*.txt"), ("All files", "*.*") ))
    CSV_Text1.delete(0.0, END)
    CSV_Text1.insert(0.0, root.CSVinputfilepath1)
    CSVinputfilename1 = root.CSVinputfilepath1[root.CSVinputfilepath1.rfind("/")+1:]
    pathCSV1 = root.CSVinputfilepath1[:root.CSVinputfilepath1.rfind("/")+1]
    CSVinput_extension1 = root.CSVinputfilepath1[root.CSVinputfilepath1.rfind("."):]
    try:
        f = open(root.CSVinputfilepath1, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_CSVfile1 = open(root.CSVinputfilepath1, 'r')
        # read the content of the file line by line 
        CSVdata_input1 = input_CSVfile1.readlines()
        #Submit input file's contents(.CSV) onto the respective text box for review by the user
        CSV_file_Text1.delete(0.0, END)
        CSV_file_Text1.insert(0.0, CSVdata_input1)
def browseMSH():
    root.MSHinputfilepath = filedialog.askopenfilename(initialdir="C://Desktop/", title = "Select a file", filetypes = ((".msh files","*.msh"), (".txt files", "*.txt"), (".csv files","*.csv"), ("All files", "*.*") ))
    MSH_Text.delete(0.0, END)
    MSH_Text.insert(0.0, root.MSHinputfilepath)
    MSHinputfilename = root.MSHinputfilepath[root.MSHinputfilepath.rfind("/")+1:root.MSHinputfilepath.rfind(".")]
    pathMSH = root.MSHinputfilepath[:root.MSHinputfilepath.rfind("/")+1]
    MSHinput_extension = root.MSHinputfilepath[root.MSHinputfilepath.rfind("."):]
    try:
        f = open(root.MSHinputfilepath, 'r')
    except:
        tkinter.messagebox.showerror(title = "Error", message = "File not found or path is incorrect")
    finally:
        input_MSHfile = open(root.MSHinputfilepath, 'r')
        # read the content of the file line by line 
        MSHdata_input = input_MSHfile.readlines()
        #Submit input file's contents(.CSV) onto the respective text box for review by the user
        MSH_file_Text.delete(0.0, END)
        MSH_file_Text.insert(0.0, MSHdata_input)
def createUpscaleMap_txt():
    CSVinputfilename1 = root.CSVinputfilepath1[root.CSVinputfilepath1.rfind("/")+1:root.CSVinputfilepath1.rfind(".")]
    pathCSV1 = root.CSVinputfilepath1[:root.CSVinputfilepath1.rfind("/")+1]
    CSVinput_extension1 = root.CSVinputfilepath1[root.CSVinputfilepath1.rfind("."):]
    savewd = pathCSV1    
    MSHinputfilename = root.MSHinputfilepath[root.MSHinputfilepath.rfind("/")+1:root.MSHinputfilepath.rfind(".")]
    pathMSH = root.MSHinputfilepath[:root.MSHinputfilepath.rfind("/")+1]
    MSHinput_extension = root.MSHinputfilepath[root.MSHinputfilepath.rfind("."):]
    try:
        # path = root.inputfilepath[:root.inputfilepath.rfind("/")+1]
        header1 = int(header_Entry1.get())
        readcols1 = [int(i) for i in readcols_Entry1.get().split(",")]
        skiprows1 = [i for i in skiprows_Entry1.get().split(",")]
        unsampled_EAV = simplest_type(UnsampledEAV_Entry.get())
        '''Global Variables'''
        MaterialGroupToUpscale = simplest_type(MaterialGroupToUpscale_Entry.get())
        rangeM = simplest_type(rangeM_Entry.get())
        rangem = simplest_type(rangem_Entry.get())
        angle = simplest_type(angle_Entry.get())
        variogram_model = simplest_type(varmdl_method.get())
        sill = simplest_type(sill_Entry.get())
        nugget = simplest_type(nugget_Entry.get())
        dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
        angOK=simplest_type(360-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html
        angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
        ratio = rangem/rangeM
        param = [sill, rangeM, nugget] #sill, range, nugget
        if "" in skiprows1:
            df = pd.read_csv(pathCSV1 + CSVinputfilename1 + CSVinput_extension1, usecols = readcols1, header = header1)
        else:
            skiprows = [simplest_type(i) for i in skiprows]
            df = pd.read_csv(pathCSV1 + CSVinputfilename1 + CSVinput_extension1, usecols = readcols1, header = header1, skiprows = skiprows1)           #  
    except:
        tkinter.messagebox.showerror(title = "Error", message = "There was an error reading the input file.\nThe most likely reason is that the file does not conform with the required format.\n\nMake sure the file has only 3 columns (x, y and aperture, respectively) and a header on row 1.")
    finally:
        if "" in skiprows1:
            df = pd.read_csv(pathCSV1 + CSVinputfilename1 + CSVinput_extension1, usecols = readcols1, header = header1)
        else:
            skiprows = [int(i) for i in skiprows]
            df = pd.read_csv(pathCSV1 + CSVinputfilename1 + CSVinput_extension1, usecols = readcols1, header = header1, skiprows = skiprows1)        
        dfcolnames1 = df.columns
        xcol1, ycol1, vcol1 = dfcolnames1
        stats=df.describe().transpose(); xmin1 = stats['min'][xcol1]; xmax1 = stats['max'][xcol1]; ymin1 = stats['min'][ycol1]; ymax1 = stats['max'][ycol1];
        vmin1 = stats['min'][vcol1]; vmax1 = stats['max'][vcol1]; # nvmin = stats['min']['N'+vcol]; nvmax = stats['max']['N'+vcol]
        df = df[[xcol1, ycol1, vcol1]] # Filter DataFrame to only the necessary columns

        MaterialGroupToUpscale = simplest_type(MaterialGroupToUpscale_Entry.get())
        rangeM = simplest_type(rangeM_Entry.get())
        rangem = simplest_type(rangem_Entry.get())
        angle = simplest_type(angle_Entry.get())
        variogram_model = simplest_type(varmdl_method.get())
        sill = simplest_type(sill_Entry.get())
        nugget = simplest_type(nugget_Entry.get())
        dic={'sill': sill, 'range': rangeM, 'nugget': nugget}
        angOK=simplest_type(-angle) # Because pykrig.ok.OrdinaryKriging takes angle values in CCW orientation, i assume from North. The documentation reads: "anisotropy_angle (float, optional) – CCW angle (in degrees) by which to rotate coordinate system in order to take into account anisotropy. Default is 0 (no rotation). Note that the coordinate system is rotated." From https://geostat-framework.readthedocs.io/projects/pykrige/en/stable/generated/pykrige.ok.OrdinaryKriging.html
        angmath = 90.0 - angle # The mathematical azimuth is measured counterclockwise from EW and not clockwise from NS as the conventional azimuth is
        ratio = rangem/rangeM
        param = [sill, rangeM, nugget] #sill, range, nugget

        dfN, dfE, ENC_df, MGs, Etype = msh_of_SplitFunction_v4(pathMSH, MSHinputfilename)
        
        if prediction_method.get() == "Element Centre Kriging":

            # =====================================================================================================================================================
            print("Kriging...")
            # =====================================================================================================================================================
            DiffMap = df.astype({xcol1: np.float16, ycol1: np.float16, vcol1: np.float16})
            data = np.array(
                [
                    DiffMap[xcol1],
                    DiffMap[ycol1],
                    DiffMap[vcol1]
                ]
                ).transpose()
            if varmdl_method.get() == "Spherical":
                OK = OrdinaryKriging(
                    data[:, 0],
                    data[:, 1],
                    data[:, 2],
                    variogram_model="spherical",
                    anisotropy_angle = angOK,
                    variogram_parameters = dic,
                    anisotropy_scaling = ratio,
                    exact_values  = True,
                    verbose=True,
                    enable_plotting=False,
                ) # scaling is the ratio between the major and minor directions' ranges
    
            ###############################################################################
            # Creates the kriged grid and the variance grid. Allows for kriging on a rectangular
            # grid of points, on a masked rectangular grid of points, or with arbitrary points.
            # (See OrdinaryKriging.__doc__ for more information.)
            # style (str) – Specifies how to treat input kriging points. Specifying ‘grid’ treats xpoints and ypoints as two arrays of x and y coordinates that define a rectangular grid.
            # Specifying ‘points’ treats xpoints and ypoints as two arrays that provide coordinate pairs at which to solve the kriging system.
            # Specifying ‘masked’ treats xpoints and ypoints as two arrays of x and y coordinates
            with Spinner():
                z, ss = OK.execute("points", ENC_df['ElementCentreX-Coordinate'], ENC_df['ElementCentreY-Coordinate'], backend = "loop")
                # z is the value, ss should stand for sigma squared which is the variance. Backend options are vectorized or loop: vectorized faster but memory intensive whereas loop is slower but also less memory-intensive.
                ###############################################################################
                                
                if type(simplest_type(unsampled_EAV)) is str and unsampled_EAV.lower() == "max":                      # unsampled Element Attribute Value
                    unsampled_EAV = np.nanmax(z) 
                elif type(simplest_type(unsampled_EAV)) is str and unsampled_EAV.lower() == "min":
                    unsampled_EAV = np.nanmin(z)
                elif type(simplest_type(unsampled_EAV)) is str and unsampled_EAV == "":
                    unsampled_EAV = nan
                else:
                    unsampled_EAV = simplest_type(unsampled_EAV) # unsampled Element Attribute Value
                
                ENC_df[vcol1] = z
                ENC_df = ENC_df.sort_values(by=['MaterialGroup', 'ElementNumber'])
            
        elif prediction_method.get() == "Element Averaging":
            # #############################       Find the nodes for each element and calculate its aperture mean      ############################# 
            dfE = dfE.sort_values(by=['MaterialGroup', 'ElementNumber'])
            dfE['eam'] = 0 # Creates a zeros column
            dfE['ElementCentreX'] = 0
            dfE['ElementCentreY'] = 0
            if "qua" in Etype.lower():
                for i, element in dfE.iterrows():
                    # e = [element['Node1'], element['Node2'], element['Node3'], element['Node4']]; # Decprecated. Hardcoded
                    e = [ element[f"Node{n+1}"] for n in range( len([n for n in dfE.columns if n.find("Node") != -1]) )]        # Grab the element's nodes' numbers
                    dfN_filt = dfN[dfN.NodeTag. isin(e)]                                                                        # Create a filtered df with only the element's nodes
                    df_filt = df[df.x.between(  dfN_filt['x'].min(), dfN_filt['x'].max() )]                                     # Create a filtered df with only the points inside the element x-coordinates
                    df_filt = df_filt[df_filt.y.between(  dfN_filt['y'].min(), dfN_filt['y'].max() )]                           # Create a filtered df with only the points inside the element y-coordinates
                    dfE.loc[i, 'ElementCentreX'] =  df_filt['x'].mean()                                                         # Calculate the element's x-coordinate centre
                    dfE.loc[i, 'ElementCentreY'] =  df_filt['y'].mean()                                                         # Calculate the element's y-coordinate centre
                    if df_filt.size==0:                                                                                         # Control in case there are no samples in within the element. We were having issues when using upscaled data because when upscaling the xy df datapoints shrink.
                        eam = eamthreshold
                    else:
                        eam = df_filt["DiffMap"].mean()                                                                         # Element Apertures' average
                        # if eam < eamthreshold:                                                                                  # Control if eam below the user-defined threshold
                        #     eam = eamthreshold
                    dfE.loc[i, 'eam'] =  eam                                                                                    # Assign element average
                    print(dfE.loc[i])
# =============================================================================
#             # Disabled for non-McDermott 2015 cases
#             dfE.loc[dfE["eam"]<=eamthreshold, "eam"] = nsmall(dfE['eam'].unique(), 2) # Replace eam column <= eamthreshold with the minimum aperture of dfE after all points are average into their respective elements. This is to manage the spread in the apertures and conductivities txt files for OGS not to complain.
# =============================================================================
            print(dfE)
            if unsampled_EAV == "max":                      # unsampled Element Attribute Value
                unsampled_EAV = dfE["eam"].max()  
            elif unsampled_EAV == "min":
                unsampled_EAV = dfE["eam"].min()
            elif unsampled_EAV == "":
                unsampled_EAV = nan
            else:
                unsampled_EAV = simplest_type(unsampled_EAV) # unsampled Element Attribute Value
            
# =============================================================================
#             # Disabled for non-McDermott 2015 cases
#             dfE.loc[dfE['MaterialGroup'] != MaterialGroupToUpscale, "eam"] = unsampled_EAV
# =============================================================================
            

        # #############################       Write upscaled attribute .txt files      ############################# 
        # Write upscaled attribute.txt
        outputfilepath = pathMSH + MSHinputfilename + "_AttributeMeshMapping_v4" + output_extension
        with open(outputfilepath, 'w') as output_file:        
            output_file.write("#MEDIUM_PROPERTIES_DISTRIBUTED\n")	
            output_file.write("$MSH_TYPE\n")
            output_file.write(" LIQUID_FLOW\n")
            output_file.write("$MMP_TYPE\n")
            output_file.write(" GEOMETRY_AREA\n")		
            output_file.write("$DIS_TYPE\n")		
            output_file.write(" ELEMENT\n")		
            output_file.write("$DATA\n")
            for mg in np.sort(dfE["MaterialGroup"].unique() ): # For each Material Group
                if type(MaterialGroupToUpscale) is str and MaterialGroupToUpscale.lower() == "all":
                    for i, element in ENC_df.iterrows():
                        if prediction_method.get() == "Element Averaging":
                            output_file.write(f"{element['ElementNumber']}\t{element['eam']}\n")
                        elif prediction_method.get() == "Element Centre Kriging":
                            output_file.write(f"{element['ElementNumber']}\t{element[vcol1]}\n")                        
                elif type(MaterialGroupToUpscale) is int and mg == MaterialGroupToUpscale:
                    for i, element in ENC_df.iterrows():
                        if element["MaterialGroup"] == MaterialGroupToUpscale and element["MaterialGroup"] != nan:
                            if prediction_method.get() == "Element Averaging":
                                output_file.write(f"{element['ElementNumber']}\t{element['eam']}\n")
                            elif prediction_method.get() == "Element Centre Kriging":
                                output_file.write(f"{element['ElementNumber']}\t{element[vcol1]}\n")
                        else:
                            output_file.write(f"{element['ElementNumber']}\t{unsampled_EAV}\n")
                elif type(MaterialGroupToUpscale) is int:
                    for i, element in ENC_df.iterrows():
                        if element["MaterialGroup"] != MaterialGroupToUpscale or element["MaterialGroup"] == nan:
                            output_file.write(f"{element['ElementNumber']}\t{unsampled_EAV}\n")
                else:
                    print(f"The provided 'Material Group To Upscale' input is of type {type(MaterialGroupToUpscale)}. The type of this input has to be int.")
            output_file.write("$STOP")
            output_file.write(f"// Averaging method: {prediction_method.get()}")
            if prediction_method.get() == "Element Centre Kriging":
                output_file.write(f"// rangeM: {rangeM}, rangem: {rangem}, angle: {angle}, variogram_model: {variogram_model}, sill: {sill}, nugget{nugget}")
                
        dfE.to_csv( pathMSH + MSHinputfilename + "_dfE.csv")
        dfN.to_csv( pathMSH + MSHinputfilename + "_dfN.csv")
        ENC_df.to_csv( pathMSH + MSHinputfilename + "_ENC_df.csv")


'''Labels'''
#Create labels
CSV_Label1 = Label(BROWSE_frame, text = "Attribute .CSV file", pady = 5)
MSH_Label = Label(BROWSE_frame, text = ".MSH file", pady = 5)
header_Label1 = Label(BROWSE_frame, text = "Header row:", pady = 5)
readcols_Label1 = Label(BROWSE_frame, text = "Columns to read (separete by comma):", pady = 5)
skiprows_Label1 = Label(BROWSE_frame, text = "Rows to skip (separete by comma):\nLeave blank if no rows to skip.", pady = 5)
csv_file_label1 = Label(CSV_frame1, text = "Attribute .CSV file content", pady = 5)
msh_file_label = Label(MSH_frame, text = ".MSH file content", pady = 5)
UnsampledEAV_Label = Label(PRED_frame, text = "Unsampled Element Attribute Value:\nOptions: max, min, na or a specific value\nLeave blank if using element average.", pady = 5)
MaterialGroupToUpscale_label = Label(PRED_frame, text = "Material Group to Upscale\nOptions: All or a specific Material Group identifier integer", pady = 5)
rangeM_label = Label(VARMODEL_frame, text = "Major continuity direction range", pady = 5)
rangem_label = Label(VARMODEL_frame, text = "Minor continuity direction range", pady = 5)
angle_label = Label(VARMODEL_frame, text = "Major contiuity direction angle (CW from North)", pady = 5)
variogram_model_label = Label(VARMODEL_frame, text = "Variogram_model (spherical)", pady = 5)
sill_label = Label(VARMODEL_frame, text = "Sill", pady = 5)
nugget_label = Label(VARMODEL_frame, text = "Nugget", pady = 5)
prediction_method_label = Label(PRED_frame, pady = 5, text = "Prediction Method")
varmdl_method_label = Label(VARMODEL_frame, pady = 5, text = "Variogram Model")

'''Text boxes and scrolls'''
#Create scrollbars
CSVscroll1 = Scrollbar(CSV_frame1)
MSHscroll = Scrollbar(MSH_frame)
#Create text boxes
CSV_Text1 = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
MSH_Text = Text(BROWSE_frame, width = 60, height=2, state=NORMAL)
header_Entry1 = Entry(BROWSE_frame, width = 10, state=NORMAL); header_Entry1.insert(-1, "0")
readcols_Entry1 = Entry(BROWSE_frame, width = 10, state=NORMAL); readcols_Entry1.insert(-1, "1,2,3")
skiprows_Entry1 = Entry(BROWSE_frame, width = 10, state=NORMAL); skiprows_Entry1.insert(-1, "")
CSV_file_Text1 = Text(CSV_frame1, height=20, state=NORMAL, yscrollcommand = CSVscroll1.set)
MSH_file_Text = Text(MSH_frame, height=20, state=NORMAL, yscrollcommand = MSHscroll.set)
UnsampledEAV_Entry = Entry(PRED_frame, width = 10, state=NORMAL); UnsampledEAV_Entry.insert(-1, "max")
MaterialGroupToUpscale_Entry = Entry(PRED_frame, width = 10, state=NORMAL); MaterialGroupToUpscale_Entry.insert(-1, "0")
rangeM_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangeM_Entry.insert(-1, "5") #rangeM_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangeM_Entry.insert(-1, "12")
rangem_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangem_Entry.insert(-1, "2") #rangem_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); rangem_Entry.insert(-1, "8")
angle_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); angle_Entry.insert(-1, "157.5") #angle_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); angle_Entry.insert(-1, "135")
variogram_model_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); variogram_model_Entry.insert(-1, "spherical")
sill_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); sill_Entry.insert(-1, "1")
nugget_Entry = Entry(VARMODEL_frame, width = 10, state=NORMAL); nugget_Entry.insert(-1, "0")
#Configure scrollbar
CSVscroll1.config(command=CSV_file_Text1.yview)
MSHscroll.config(command=MSH_file_Text.yview)

'''Buttons'''
BrowseCSVButt1 = Button(BROWSE_frame, text="Browse Attribute .CSV file", fg="black", font=("Ariel", 9, "bold"), command=browseCSV1)
BrowseMSHButt= Button(BROWSE_frame, text="Browse .MSH file", fg="black", font=("Ariel", 9, "bold"), command=browseMSH)
FormatButt= Button(BROWSE_frame, text="Create attribute to Gmsh/GINA .msh element Mapping .txt file", fg="black", font=("Ariel", 9, "bold"), command=createUpscaleMap_txt)
"Radio Buttons"
avg_RButt = Radiobutton(PRED_frame, text = "Element Averaging", variable = prediction_method, value = "Element Averaging", command = prediction_method_selection)
krig_RButt = Radiobutton(PRED_frame, text = "Element Centre Kriging", variable = prediction_method, value = "Element Centre Kriging", command = prediction_method_selection)
sph_varmdl_RButt = Radiobutton(VARMODEL_frame, text = "Spherical", variable = varmdl_method, value = "Spherical", command = varmdl_method_selection)

'''Allocate widgets'''
# BROWSE_frame
BrowseCSVButt1.grid(row = 0, column = 2, sticky = W)
BrowseMSHButt.grid(row = 0, column = 5, sticky = W)
header_Label1.grid(row = 2, column = 0, sticky = W)
header_Entry1.grid(row = 2, column = 1, sticky = W)
readcols_Label1.grid(row = 3, column = 0, sticky = W)
readcols_Entry1.grid(row = 3, column = 1, sticky = W)
skiprows_Label1.grid(row = 4, column = 0, sticky = W)
skiprows_Entry1.grid(row = 4, column = 1, sticky = W)
FormatButt.grid(row = 0, column = 6, sticky = W)
# MSH_frame
MSH_Label.grid(row = 0, column = 3, sticky = W)
MSH_Text.grid(row = 0, column = 4, sticky = W)
msh_file_label.grid(row = 1, column = 1, sticky = W)
MSH_file_Text.grid(row = 2, column = 1)
MSHscroll.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
#CSV_frame
CSV_Label1.grid(row = 0, column = 0, sticky = W)
CSV_Text1.grid(row = 0, column = 1, sticky = W)
csv_file_label1.grid(row = 1, column = 1, sticky = W)
CSV_file_Text1.grid(row = 2, column = 1)
CSVscroll1.grid(row = 1, column = 3, rowspan = 3, sticky = NS)
# PRED_frame
MaterialGroupToUpscale_label.grid(row = 0, column = 0, sticky = W)
MaterialGroupToUpscale_Entry.grid(row = 0, column = 1, sticky = W)
prediction_method_label.grid(row = 1, column = 0, sticky = W)
avg_RButt.grid(row = 1, column = 1, sticky = W)
krig_RButt.grid(row = 2, column = 1, sticky = W)
UnsampledEAV_Label.grid(row = 4, column = 0, sticky = W)
UnsampledEAV_Entry.grid(row = 4, column = 1, sticky = W)
#VARMODEL_frame
varmdl_method_label.grid(row = 0, column = 0, sticky = W)
sph_varmdl_RButt.grid(row = 0, column = 1, sticky = W)
rangeM_label.grid(row = 1, column = 0, sticky = W)
rangem_label.grid(row = 2, column = 0, sticky = W)
angle_label.grid(row = 3, column = 0, sticky = W)
sill_label.grid(row = 4, column = 0, sticky = W)
nugget_label.grid(row = 5, column = 0, sticky = W)
rangeM_Entry.grid(row = 1, column = 1, sticky = W)
rangem_Entry.grid(row = 2, column = 1, sticky = W)
angle_Entry.grid(row = 3, column = 1, sticky = W)
sill_Entry.grid(row = 4, column = 1, sticky = W)
nugget_Entry.grid(row = 5, column = 1, sticky = W)

'''Menu bar'''
menubar = Menu(root)
root.config(menu = menubar)

file_menu = Menu(menubar)
menubar.add_cascade(label = "File", menu = file_menu)
file_menu.add_command(label = "Create attribute to Gmsh/GINA .msh element Mapping .txt file", command = createUpscaleMap_txt)
file_menu.add_separator()
file_menu.add_command(label="Exit", command=root.quit)

options_menu = Menu(menubar)
menubar.add_cascade(label = "Options", menu = options_menu)
options_menu.add_command(label = "Description", command = description)

'''Prompt on open'''
#Uncomment below if you wish the description() function (info box explaining how the app works) to be prompted at app opening
#description()

'''Add Window title, geometry and create Window's main loop'''
root.title("Tkinter - CSV attribute to Gmsh/OGS-GINA .msh element Mapping v4")
root.geometry("2000x800")
root.mainloop()