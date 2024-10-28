# QM vs MM geometry parameter dashboard

This code calculates the geometric paramenters (e.g. bond lenght, angles, dihedral angles) for a set of molecules and compares the values between the QM-optimized molecule and MM-optimized molecule.

## Python environment
The conda environment used to develop the code is included as `env.yaml`. This is just a minimal environment needed to work with this code, if you need other packages, please add them.

## Calculating the geometric parameters
There are two ways you can analyze molecular geometries: either using a database produced by our [internal YAMMBS benchmarking code](https://github.com/openforcefield/yammbs), or by reading geometries from SDF files.

### Reading from YAMMBS database
If reading in a database generated by YAMMBS, you need to use the `--db` flag to read in the database, which can be located anywhere. 
The `--ff_yammbs` flag can be used to identify which force field's results you want to use within the database, however it will default to using the first force field stored.
This is only used for fetching records, and the actual force field does not need to be present.

### Reading from SDF files
First you must have two directories of molecular structures, one with SDF files for QM-optimized geometries and one with SDF files for MM-optimized geometries, specified by `--qm_dir`, `--mm_dir`, as described below.
They can be in an arbitrary location, passed with the `--dir` flag.
Within each directory, molecules should be named according to the convention mol-i-conf-j.sdf, where i is the molecule number (can be anything) and j is the two-digit conformer number within that molecule (e.g. 00, 01, 02...).
If using the default option to only analyze one conformer per molecule (`--conformers False`), the code expects each molecule to have `mol-i-conf-00.sdf`, or that molecule will be skipped. 
SDF files must have a field called `Record QCArchive` with the QCArchive ID (or some other ID) if you want to use the `problem_file` option to filter out problematic QCArchive IDs.

### Calculation options
If using a YAMMBS database, the following options are needed:

`--db`: YAMMBS SQlite database to read from (if using)

`--ff_yammbs`: Force field used to calculate geometries with yammbs (if reading from database). If not provided, it will use the first stored force field in the database.


If reading from SDF files, the following options are needed:

`--mm_dir`: Directory with MM optimization SDF files (if not reading from database)

`--qm_dir`: Directory with QM optimization SDF files (if not reading from database)

`--dir`: Directory where QM and MM directories are located. Defaults to the working directory.


For both cases, the following options are needed:

`--ff_file`: Force field to group parameters by. If not a released force field, must be the full path to the file. (Default: openff-2.1.0.offxml)

`--label`: Label to save data files with. Data will be stored in a directory with this name, with files of the name `angles_qmv[label].json`, `bonds_qmv[label].json, etc.


And optionally:

`--conformers`: Whether to use all conformers. If False, just use the first conformer for each molecule. Note that setting this to True can lead to very large files and long run times.

`--filter`: (Optional) SMARTS pattern to filter for. Must have at least one tagged atom.

`--problem_file`: (Optional) File(s) listing QCAIDs of conformers to exclude from analysis.

### Running the code
You can then calculate the parameters for the two datasets using:

`python calculate_params.py [flags]`

This will create four data files, one for bond lengths, one for angles, one for "Proper torsion" dihedral angles, and one for "Improper torsion" dihedral angles.

Examples using both YAMMBS databases and SDF files with a variety of options can be found in the `example` directory.

## Visualizing the parameters
Once you have these data files, you can visualize them by running the dashboard:

`python dashboard.py [data directory]`

This will create a local server. To use the dashboard, go to the address printed (e.g. `Dash is running on [address]`)

Once the dashboard is running, you can toggle between the different data types and the different force field parameters. 
If you click on the data points, it will show the structure of that molecule, with the parameter in question highlighted.

![Screen Shot 2024-01-04 at 6 00 21 PM](https://github.com/amcisaac/geom_dash/assets/29759281/91203eb0-1446-4ccc-9ec3-ffe3ca7d80eb)

Additionally, you can filter the data based on a SMIRKs pattern, to identify trends in the data. In the box that says "Enter a tagged SMIRKs", enter a pattern, where the appropriate number of atoms are tagged (e.g. for a bond, 2 atoms tagged, for an angle, 3 atoms tagged). Molecules that match the SMIRKs pattern will be highlighted in red.

![Screen Shot 2024-01-04 at 5 58 21 PM](https://github.com/amcisaac/geom_dash/assets/29759281/abb8d1e9-f172-440c-82d2-5d4e2dc8f2cc)


## Note on the code

Some of the dashboard code was borrowed from the OpenFF docs or from Brent's Lipoma repo [here](https://github.com/ntBre/lipoma). These instances are noted in the comments of the code.
