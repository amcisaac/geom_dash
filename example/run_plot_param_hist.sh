#!/bin/bash

# Plot one benchmark
python ../plot_param_hist.py --data_dir sage220-fromsqlite-noconfs-noprobs --labels 'One conf, all recs' --ff_file 'openff_unconstrained-2.2.0.offxml' --bonds b1 --bonds b2 --angles a1 --angles a2 --propers t1 --propers t2 --impropers i1 --impropers i2 --save_dir param_plots_sage220-fromsqlite-noconfs-noprobs

# Plot two benchmarks (usually this would be two different FFs)
python ../plot_param_hist.py --data_dir sage220-fromsqlite-noconfs-noprobs --labels 'One conf, all recs' --data_dir sage220-fromsqlite-confs-noprobs --labels 'All confs, all recs' --ff_file 'openff_unconstrained-2.2.0.offxml' --bonds b1 --bonds b2 --angles a1 --angles a2 --propers t1 --propers t2 --impropers i1 --impropers i2 --save_dir param_plots_sage220-fromsqlite-noconfs-confs-noprobs
