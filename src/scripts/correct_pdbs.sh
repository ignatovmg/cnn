#!/bin/bash

cd sandbox/aligned_pdb

mkdir ../tmp

for pdb in `find ./ -name "*.pdb"` 
do
	cp $pdb ../tmp
	../../src/scripts/prepare_pdb.sh ../tmp/$pdb
	
	echo $pdb
	name=$(python -c "print('${pdb}'[2:-4])")
	mv ../tmp/${name}_nmin.pdb ../corrected_pdb/$(python -c "print('${name}'[:4])").pdb
	
	rm -f ../tmp/*
done

cd -


