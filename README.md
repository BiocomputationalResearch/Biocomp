## Project Name: Structure Deviation

### Version: 1.1

### Purpose of the Project:

To compare a coordinating site (comparable) with the standard site of the same coordination geometry as that of the comparable one on the structure of bio-molecules at a time in order to produce minimum mean square errors due to their structural deviations in sides and angles. Before comparison, the original coordinating atoms of the standard and the comparable sites are transformed by applying some mathematical transformation to make them structurally compatible.



### Algorithms Used:

Four algorithms have been used to compare pair of geometries for standard as well as comparable  viz. **1) trigonal planar/t-shaped, 2) square/tetragonal planar, 3) tetrahedral,** and **4)** to compare pair of geometries **square/tetragonal pyramidal, trigonal bipyramidal** and **octahedral/tetragonal bipyramidal** .



### Provisions for Metal Sites to be Compared:

1. For **Trigonal Planar/T-shaped Sites:** Na, Ca, K, Mg, Mn, Zn, Fe, Cu, Cd, Co, As, Pb, Hg, Ni

2. For **Square/Tetragonal Planar Sites:** Mg, Fe, Ni, Cd

3. For **Tetrahedral Sites:** Na, Ca, K, Mg, Mn, Zn, Fe, Cu, Cd, Co, As, Pb, Hg, Ni, Cs

4. For **Square/Tetragonal Pyramidal Sites:** Na, Ca, Mg, Zn, Cd, Ni

5. For **Trigonal Bipyramidal Sites:** Na, Ca, Mg, Mn, Zn, Fe, Cd, Co, Hg, Ni

6. For **Octahedral/Tetragonal Bipyramidal Sites:** Na, Ca, Mg, Mn, Zn, Fe, Cd, Co, Ni



### Tasks to be performed:

1. Compares batch of input structures with the respective standard structure after making them structurally compatible.

2. Produces results of comparison for each pair of geometries in report form. 

3. Generates results of distance and angular deviations as individual report file for each metal site present in an input structure

4. Generates a summary report file for each type of metal present in all comparable input structures.

5. Generates diagrams to visualise the concept used in the algorithms.

6. Extracts coordinating atoms for each metal site in .pdb format to be viewed in Rasmol for both of the standard and the comparable sites.

7. Computes results of both types of deviations obtained with original coordinating atoms.

   

### Input Format: .pdb or .cif

### Output Format: 

* .mbsi - for report files
* .pdf - for illustrations
* .stdn - coordinating atoms with bond information for standard structure in .pdb format - for original and transformed atoms
* .cmp - coordinating atoms with bond information for comparable structure in .pdb format - for original and transformed atoms
* .comp - Deviations output obtained with original atoms



### Programming language Used: C++

### Libraries Used: Two C++ Libraries Openbabel and Gemmi

* Openbabel is required to convert .pdb/.cif format data to .mol2 format prior to any computation in order to obtain bond information between any two atoms.

* Gemmi is required to read .cif file.

  

### Compilation & Execution:

* Prior to compilation, **Openbabel and Gemmi are needed to be installed**.

* After compilation, while executing the program **except the input filenames along with extension** , normally **no other options are required**.

* It has been noticed that with input file whose size is upto 3 MB, number of bonds is within 20000 and number of different kind of metals present is 5-6 each, this Project is working satisfactorily. But with large file having huge numbers of bonds and more numbers of metals in it, sometimes the program crashes. So to save the program from being crashed, while executing the program an option (-M) is used to set a threshold value for number of metals which may vary from system to system. Normally, a default value 5 is used as this threshold value, when no option is used with the command of execution. 

* The commands to compile & execute the program:

  cmake ..

  make

  ./Structure_Deviation 2a2a.pdb [2a2a.cif] - with no option

  ./Structure_Deviation -M 8 2a2a.pdb [2a2a.cif] - with option

