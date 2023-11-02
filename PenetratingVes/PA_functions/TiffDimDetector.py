#!/usr/bin/env python
# coding: utf-8

# In[11]:


# CREATED: 7-DEC-2022
# LAST EDIT: 8-DEC-2022
# AUTHOR: DUANE RINEHART, MBA (drinehart@ucsd.edu)
# Jacob Duckworth (jaduckwo@ucsd.edu)

# Read tif file with python-bioformats

from pathlib import Path, PurePath, PureWindowsPath
from tifffile import imread,TiffFile
import tifffile
import sys
#################################################################
input_path = Path(sys.argv[1])
input_path = str(input_path).replace("**"," ")
#################################################################
# print(f"argumentcaptured:{input_path}")
print(input_path)

tif = TiffFile(input_path) #Commented out

print(len(tif.pages))
series = tif.series[0];
print(series.shape)

numims = len(tif.pages)
