# 3-D-Deformation-Inversion
<img width="416" alt="image" src="https://github.com/arvin1367/3-D-Deformation-Inversion/assets/49364261/5d3bd8bd-1613-4d0a-ad6c-36321716a209">

What is this repository for?

A MATLAB toolbox for automatically calculating SAR-derived 3-D deformation maps of glacier, landslide, and land subsidence

Please cite the software as: Oriani F., Treble P. C., Baker A., Mariethoz G., WlCount: Geological lamination detection and counting using an image analysis approach, Computers & Geoscience, https://doi.org/10.1016/j.cageo.2022.105037

How do I get set up?
See the tutorial video: https://odysee.com/@Fabio:d7/wlcount_tutorial:3?r=CrZbidUqeeBLeyCokTWSSM6bemSKux1p or follow the instructions: Copy the repository folder ( Downloads from the menu on the left of this page) in Windows 64bit. In alternative use any other machine with python3 installed. If using the non-compiled python script, satisfy the dependencies listed.

Usage
1) Choose one of the following:

For Win64 users: go in the dist directory and double click on wlcount_gui.exe or, from a terminal, cd to wlcount/dist and type wlcount_alpha_gui.exe.
For other users with python3: from a terminal launch python wlcount_gui.py.
2) After some seconds a dialog window opens to select the image to be analysed. For testing purposes choose test_image.tif from the wlcount_alpha folder.

3) A dialog window appears asking if you want to use a previous allignment: if it is the first time you perform the analysis in this location, choose no. Otherwise choose yes to reuse the image allignment files previously created and save computational time (see OUTPUT FILES below).

4) When the image appears, draw a west-east section in the wavelet image to perform the count. To do this, draw at least 2 pilot points from left to right in the wavelet image (bottom graph) using the mouse pointer:

left click = draw a new point
right click = cancel the previous point
ctrl+r = cancel and restart the current section
enter button = automatically draw the section among the points and perform the count. The count appears in the central image legend.
5) repeat point 4) to make other counts. ctrl+e to end the analysis or close the image to abort the count.

6) The program ends displaying the counts in the terminal and creating the output files in the image folder.

OUTPUT FILES (in the same image folder):

<image_name>.pdf --> display the image of the counts
<image_name>.npy --> image allignment file (reusable)
<image_name>_pic.npy --> image allignment file (reusable)
<image_name>_counts.dat --> dat file containing the 0/1 time series of all counts performed. They have the same length of the x axis of the image: 1 corresponds to the center of a lamina counted and 0 everywhere else.
