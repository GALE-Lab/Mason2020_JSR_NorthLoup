# Mason2020_JSR_NorthLoup
Subaqueous dune field pattern evolution and interactions: North Loup River, Nebraska, USA 
Mason et al. in review. 
For more information see the full publication or contact the corresponding author at daym@epss.ucla.edu

This repository includes data used in the publication with the above title, currently in review with the Journal of Sedimentary Research. 

Contents of this repository include: 
1) Images with associated ancillary files for georeferencing those images 
2) Shapefiles of dune crestlines and spurs traced and interpreted from images 
3) An xls table with x and y point data associated with images and dune crests 
4) A Matlab script used to conduct the associated data analysis and quantification of deformation 

1) Images: 
The research conducted in this work used time-lapse images collected by a camera suspended over the North Loup river, taking pictures every minute for several hourse. The image files included all begin with DSC_ the subsequent number in the filename equates to time in minutes, with the study beginning at t0 = 5869 or DSC_5869.JPG. JPG files with additional extensions can be used to georeference the images in GIS platforms. 

2) Shapefiles: 
Dune and spur position was mapped on each of the images and is recorded as a shapefile compatible with GIS platforms. Files begin with day3morning_ followed by the time stamp that matches the image in which the mapping was conducted. 1_200_scale refers to the scale at which tracing was conducted. 

3) XLS Table 
The XLS table of point data includes one sheet per frame and records in each row an ID number for the crest on which the point was mapped and the x and y position of that point in UTM coordinates. 

4) MATLAB Script 
The data analysis included in the main text of this publication calculated deformation of mapped crestlines between time steps. The Matlab script in this repository was used to make these measurements. See comments in the code itself for more details. 

For more information, contact Mackenzie Day at daym@epss.ucla.edu
