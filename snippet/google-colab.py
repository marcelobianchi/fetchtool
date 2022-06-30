#
# This snippet is a tricky to get obspy up and running on google 
# colab and on the sequence install fetchtool on it.
#
# last tested -- 2022-06-30
#
##### 

##
# This is a tricky to get cartopy up and running into the Colab Notebook
# https://colab.research.google.com/github/astg606/py_materials/blob/master/visualization/introduction_cartopy.ipynb
##

!apt-get install libproj-dev proj-data proj-bin libgeos-dev python3-cartopy python-cartopy cython
!pip install cartopy==0.18.0
!pip uninstall -y shapely
!pip install shapely --no-binary shapely

##
# Here is the obspy
##

!pip install obspy

## -- break a cell here --

#
# After obspy is installed you need to restart the environment, 
# google colab will show you a button to click on it.
#

##
# Here we install fetchtool
##

!pip install git+https://github.com/marcelobianchi/fetchtool.git
