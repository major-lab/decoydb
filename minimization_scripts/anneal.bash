#!/bin/bash
# $1 gzipped pdb file to anneal

pdb_file=$1
pdb_prefix=`echo $pdb_file | cut --delimiter='.' -f 1`

/u/major/bin/Tinker/addions.exe $pdb_prefix.pdb
/u/major/bin/Tinker/get_sched.exe $pdb_prefix.pdb
/u/major/bin/Tinker/prepchains.exe $pdb_prefix.pdb
/u/major/bin/Tinker/pdbxyz $pdb_prefix.pdb ALL /u/major/bin/Tinker/amber99
/u/major/bin/Tinker/anneal $pdb_prefix.xyz -k /u/major/bin/Tinker/gbsa.key 100 0 0 500 L 1.0 0.5 0
/u/major/bin/Tinker/xyzpdb $pdb_prefix.001 /u/major/bin/Tinker/amber99
rename pdb_2 pdb.min *.pdb_2
/soft/bioinfo/linux/mcrms-dev/bin/MC-RMSD -n $pdb_prefix.pdb $pdb_prefix.pdb.min
rename pdb.min pdb *.pdb.min
/u/major/bin/Tinker/three2one.exe $pdb_prefix.pdb
/u/major/bin/Tinker/set_sched.exe $pdb_prefix.pdb