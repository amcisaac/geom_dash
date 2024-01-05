# QM vs MM geometry parameter dashboard

This code calculates the geometric paramenters (e.g. bond lenght, angles, dihedral angles) for a set of molecules and compares the values between the QM-optimized molecule and MM-optimized molecule.

## Python environment
The conda environment used to develop the code is included as `env.yaml`. This was taken from [here](https://github.com/ntBre/lipoma) and probably not all of the packages are needed.

## Calculating the geometric parameters
To use, first you must have two directories of molecular structures, one with SDF files for QM-optimized geometries and one with SDF files for MM-optimized geometries. 
The files within the folders should have the same name for the same molecule, e.g. for `Molecule1` there should be `qm_structures/Molecule1.sdf` and `mm_structures/Molecule1.SDF`.
The path of the files and the force field is currently hard-coded.

You can then calculate the parameters for the two datasets using:

`python calculate_params.py`

This will create four data files, one for bond lengths, one for angles, one for "Proper torsion" dihedral angles, and one for "Improper torsion" dihedral angles.

## Visualizing the parameters
Once you have these datafiles, you can visualize them by running the dashboard:

`python dashboard.py`

This will create a local server. To use the dashboard, go to the address printed (e.g. `Dash is running on [address]`)

Once the dashboard is running, you can toggle between the different datatypes and the different force field parameters. 
If you click on the data points, it will show the structure of that molecule, with the parameter in question highlighted.

![Screen Shot 2024-01-04 at 6 00 21 PM](https://github.com/amcisaac/geom_dash/assets/29759281/91203eb0-1446-4ccc-9ec3-ffe3ca7d80eb)

Additionally, you can filter the data based on a SMIRKs pattern, to identify trends in the data. In the box that says "Enter a tagged SMIRKs", enter a pattern, where the appropriate number of atoms are tagged (e.g. for a bond, 2 atoms tagged, for an angle, 3 atoms tagged). Molecules that match the SMIRKs pattern will be highlighted in red.

![Screen Shot 2024-01-04 at 5 58 21 PM](https://github.com/amcisaac/geom_dash/assets/29759281/abb8d1e9-f172-440c-82d2-5d4e2dc8f2cc)


## Note on the code

Some of this code was borrowed from the OpenFF docs or from Brent's Lipoma repo [here](https://github.com/ntBre/lipoma). These instances are noted in the comments of the code.
