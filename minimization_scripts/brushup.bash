#!/bin/bash
# $1 gzipped pdb file to brushup

pdb_file=$1
pdb_prefix=`echo $pdb_file | sed 's/^\(.\{1,\}\)\.pdb$/\1/'`

/u/major/bin/Tinker/addions.exe $pdb_prefix.pdb
/u/major/bin/Tinker/get_sched.exe $pdb_prefix.pdb
/u/major/bin/Tinker/prepchains.exe $pdb_prefix.pdb
/u/major/bin/Tinker/pdbxyz $pdb_prefix.pdb ALL /u/major/bin/Tinker/amber99
/u/major/bin/Tinker/prepkey.exe $pdb_prefix.xyz 1000 amber99 | sed 's/parameters .\{1,\}amber99/parameters \/u\/major\/bin\/Tinker\/amber99/' > $pdb_prefix.min.key
/u/major/bin/Tinker/minimize $pdb_prefix.xyz -k $pdb_prefix.min.key 0.1
/u/major/bin/Tinker/xyzpdb $pdb_prefix.xyz_2 /u/major/bin/Tinker/amber99
mv $pdb_prefix.pdb_2 $pdb_prefix.pdb.min
/soft/bioinfo/linux/mcrms-dev/bin/MC-RMSD -n $pdb_prefix.pdb $pdb_prefix.pdb.min
mv $pdb_prefix.pdb.min $pdb_prefix.pdb
/u/major/bin/Tinker/three2one.exe $pdb_prefix.pdb
/u/major/bin/Tinker/set_sched.exe $pdb_prefix.pdb
rm $pdb_prefix.seq*
rm $pdb_prefix.xyz*
rm $pdb_prefix.min.key
