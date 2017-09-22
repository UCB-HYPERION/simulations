globalSignalVis
==============================================================

This code aims to simulate the visibilities observed by the HYPERION instrument 
-- i.e. a closely packed low-frequency interferometer with absorber between 
elements. 

System Requirements
-----------------

In order to run this simulation, the following Python packages need to be 
installed:

1) Aipy
2) UVTools
3) NumPy
4) MatPlotLib

Additionally, the file absorber.py will need to be in the same directory. 

Finally, a calfile will need to be specified and in your Python path.

Usage
-----------------

To run the simulation, run the following command:

    python globalSignalVis.py --calfile name_of_calfile --simdir 
    path/to/where/simulation/data/will/be/saved 

Future Plans
-----------------

There are many improvements that we aim to make to this simulation in the 
coming weeks and months.

1) Add multi-frequency capabilities for the sky and absorber.
2) Add realistic beams for the HYPERION fat dipoles.
3) Smooth out transition between absorber and sky so it isn't a harsh cutoff.
4) Import realistic global sky models.
