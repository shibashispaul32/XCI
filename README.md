# XCI
The wave_pinning_1d directory contains more than one .in files each corresponding to different scenarios. To run the codes successfully one needs to change the name of the input file to input.in and run the code to get the simulation result for that scenario. 

The fortran codes require a fortran compiler to be installed in the system. We prefer using intel fortran compiler which is available free of cost (guide: https://gist.github.com/SomajitDey/aeb6eb4c8083185e06800e1ece4be1bd).

In addition to the Fortran codes, there are Jupyter notebooks that use Python to analyze the simulated and experimental data. To run these, you must install Python with the appropriate packages, as well as Gnuplot, on your system. For the convenience of the reader, we are providing the conda environment file which we have used to run the jupyter notebooks (find environment.yml file).

To execute the code within the Jupyter notebooks, you must first download the data folder from the following link: https://drive.google.com/file/d/16dCaZlECGkhGIXuHNL44OQbtqFrRsWOa/view?usp=sharing
After downloading, extract the contents of the .zip file to the root directory of the project. The Jupyter notebook codes can only be executed successfully after this step.

The environment.yml file can be used to to recreate the conda environment that was used in this project to execute codes.
