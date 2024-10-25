#!/bin/bash

# Only first conformer used, nothing filtered out
python ../calculate_params.py --qm_dir OpenFF_Gen2_Coverage_sage220_QM --mm_dir OpenFF_Gen2_Coverage_sage220_MM --conformers False --ff_file openff_unconstrained-2.2.0.offxml --label sage220-fromfile-noconfs-noprobs

# Only first conformer used, filter out random ids 
python ../calculate_params.py --qm_dir OpenFF_Gen2_Coverage_sage220_QM --mm_dir OpenFF_Gen2_Coverage_sage220_MM --conformers False --ff_file openff_unconstrained-2.2.0.offxml --problem_file problem_ids/random_qcaids.txt --label sage220-fromfile-noconfs-probs

# All conformers used, nothing filtered out
python ../calculate_params.py --qm_dir OpenFF_Gen2_Coverage_sage220_QM --mm_dir OpenFF_Gen2_Coverage_sage220_MM --conformers True --ff_file openff_unconstrained-2.2.0.offxml --label sage220-fromfile-confs-noprobs
 
# All conformers used, filter out random ids 
python ../calculate_params.py --qm_dir OpenFF_Gen2_Coverage_sage220_QM --mm_dir OpenFF_Gen2_Coverage_sage220_MM --conformers True --ff_file openff_unconstrained-2.2.0.offxml --problem_file problem_ids/random_qcaids.txt --label sage220-fromfile-confs-probs 
