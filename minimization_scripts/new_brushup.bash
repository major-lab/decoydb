#!/bin/bash
# $1 gzipped pdb file to brushup

pdb_file=$1
pdb_prefix=`echo $pdb_file | sed 's/^\(.\{1,\}\)\.pdb$/\1/'`

/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/addions.exe $pdb_prefix.pdb
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/get_sched.exe $pdb_prefix.pdb
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/prepchains.exe $pdb_prefix.pdb
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/pdbxyz $pdb_prefix.pdb ALL /u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/amber99
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/prepkey.exe $pdb_prefix.xyz 1000 amber99 | sed 's/parameters .\{1,\}amber99/parameters \/u\/major\/bin\/Tinker\/amber99/' > $pdb_prefix.min.key
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/minimize $pdb_prefix.xyz -k $pdb_prefix.min.key 0.1
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/xyzpdb $pdb_prefix.xyz_2 /u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/amber99
mv $pdb_prefix.pdb_2 $pdb_prefix.pdb.min
/soft/bioinfo/linux/mcrms-dev/bin/MC-RMSD -n $pdb_prefix.pdb $pdb_prefix.pdb.min
mv $pdb_prefix.pdb.min $pdb_prefix.pdb
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/three2one.exe $pdb_prefix.pdb
/u/leongs/MCSYM-soft/cgi-bin/Tinker_copy/Tinker/set_sched.exe $pdb_prefix.pdb
rm $pdb_prefix.seq*
rm $pdb_prefix.xyz*
rm $pdb_prefix.min.key
