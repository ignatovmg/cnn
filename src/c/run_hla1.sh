#!/bin/bash

gcc spring_min_hla1.c -O3 -flto -lmol2 -lm -L ~/lib -I ~/include -o spring_min_hla1

#./spring_min_hla1 spring_min_hla1_test/mutated0_nmin.pdb spring_min_hla1_test/mutated0_nmin.psf charmm_param.prm charmm_param.rtf spring_min_hla1_test/fixed0_min1.pdb spring_min_hla1_test/fixed0_min2.pdb spring_min_hla1_test/sfile spring_min_hla1_test/minimized.pdb

#./spring_min_hla1 mutated3_nmin.pdb mutated3_nmin.psf charmm_param.prm charmm_param.rtf fixed_mhc.pdb fixed_mhc.pdb sfile minimized.pdb
