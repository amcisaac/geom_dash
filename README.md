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

<img width="1246" alt="Screen Shot 2024-01-04 at 4 58 24 PM (2)" src="https://github.com/amcisaac/geom_dash/assets/29759281/53662a75-2554-48d1-913b-4e47e7b9dcfc">

## Note on the code

Some of this code was borrowed from the OpenFF docs or from Brent's Lipoma repo [here](https://github.com/ntBre/lipoma). These instances are noted in the comments of the code.
