# Octagonal Lattice Generator

This repo contains the code for generating lattice-points for $C_8$ quasicrystals.
These crystals can be tuned to host any different topological charges.


## Compiling the program

The program is compiled with ```make``` and runs automatically.
Program requires the [SFML graphical library](https://www.sfml-dev.org/) to be installed on the system
The initial geometry needs to be defined in the main.cpp -file.


## Running the program

The program is run with ```./main``` and it opens a graphical window.
![def](screenshots/default_screen.png)
Depending on the geometry, different operations can be done.
If either Ammann-Beenker or Segment are used, the system can be subdivided by pressing ```S```.
When system is subdivided into regions of sufficient size, the field can be optimized for irreps A,B,E,F,G by pressing the corresponding key.
After subdivision the points can be saved with ```W```.


### Controls
- ```S``` - Subdivides the geometry
- ```A``` - Optimizes geometry for irrep A
- ```B``` - Optimizes geometry for irrep B
- ```E``` - Optimizes geometry for irrep E
- ```F``` - Optimizes geometry for irrep F
- ```G``` - Optimizes geometry for irrep G
- ```W``` - Saves the points into *irrep*_points.csv


