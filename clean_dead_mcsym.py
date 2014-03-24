import os

decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/decoy"

pbs_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/pbs"

for elem in xrange(5, 13):
    pbs_subdir = os.path.join(pbs_dir, str(elem), "pbs")
    for pbs in sorted(os.listdir(pbs_subdir)):
        acc = pbs.replace(".pbs", "")
        print "rm -rf " + os.path.join(decoy_dir, acc)