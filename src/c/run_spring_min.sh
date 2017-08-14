#!/bin/bash

gcc spring_min.c -lmol2 -lm -o spring_min
./spring_min mutated3_nmin.pdb mutated3_nmin.psf charmm_param.prm charmm_param.rtf fixed_mhc.pdb sfile minimized.pdb

