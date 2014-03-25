import os

processed_decoy_dir_list = ["/u/leongs/reproduction_projet_naim/rel20/3D/processed/decoy",
                            "/u/leongs/reproduction_projet_naim/rel20/3D/processed_23_jan/decoy"]

backup_dir = "/export/home/leongs/processed_decoy_23_jan"

list_len = []

dict_mir = dict()

for decoy_dir in processed_decoy_dir_list:
    for acc in os.listdir(decoy_dir):
        out_dir_cont = os.listdir(os.path.join(decoy_dir, acc, 'out'))
        if not acc in dict_mir:
            dict_mir[acc] = []
        dict_mir[acc].extend([dict(path=os.path.join(decoy_dir, acc, 'out', pdb),
                                   pdb=pdb) for pdb in out_dir_cont])

for acc, list_cont in dict_mir.iteritems():
    if len(list_cont) < 1:
        print "mkdir {new}".format(new=os.path.join(backup_dir, acc))
        print "mkdir {out}".format(out=os.path.join(backup_dir, acc, "out"))
#         print "\n".join(["mv {source} {dest}".format(source=pdb_dict["path"],
#                                                      dest=os.path.join(backup_dir, acc, "out", pdb_dict["pdb"])) for pdb_dict in list_cont])
        for decoy_dir in processed_decoy_dir_list:
            if acc in os.listdir(decoy_dir):
                print "rm -rf {f}".format(f=os.path.join(decoy_dir, acc))

# for decoy in sorted(list_decoys):
#     out_dir = os.path.join(processed_decoy_dir, decoy, "out")
#     if len(os.listdir(out_dir)) < 1:
#         print "mv {source} {dest}".format(source=os.path.join(os.path.join(processed_decoy_dir, decoy)),
#                                           dest=os.path.join(backup_dir, decoy))

#     list_len.append(len(os.listdir(out_dir)))
# 
# print float(sum(list_len)) / float(len(list_len))
# print sorted(list_len)[len(list_len)//2]