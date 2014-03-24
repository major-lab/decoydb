import os

print "cd /u/leongs/reproduction_projet_naim/rel20/3D/failed_processed"

for elem in os.listdir("/u/leongs/reproduction_projet_naim/rel20/3D/failed_processed"):
    if elem.endswith(".pdb"):
        print "gzip " + os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/failed_processed", elem)