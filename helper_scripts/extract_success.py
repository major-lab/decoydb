import os
import cPickle
import shutil

success_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/success"
digested_data_pk = "/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk"
flashfold_dir = "/u/leongs/reproduction_projet_naim/rel20/2D/filtered"

dict_digested_data = dict()
with open(digested_data_pk, 'rb') as dd:
    dict_digested_data = dict((hairpin_dict["accession"], hairpin_dict) for hairpin_dict in cPickle.load(dd))

for proc in ["/u/leongs/reproduction_projet_naim/rel20/3D/processed_11_feb/decoy",
             "/u/leongs/reproduction_projet_naim/rel20/3D/processed_17_feb/decoy",
             "/u/leongs/reproduction_projet_naim/rel20/3D/processed_19_feb/decoy",
             "/u/leongs/reproduction_projet_naim/rel20/3D/processed_21_feb/decoy",
             "/u/leongs/reproduction_projet_naim/rel20/3D/processed_24_feb/decoy"]:
    for acc in os.listdir(proc):
        mcfold_output = os.path.join(flashfold_dir, acc)
        with open(mcfold_output, 'rb') as mcfold_o:
            list_struct = [elem.strip() for elem in mcfold_o]

        out_dir = os.path.join(proc, acc, "out")
        for out_pdb in os.listdir(out_dir):
            index = out_pdb.split("-")[0].split("_")[1]
            struct = list_struct[int(index)].split()[0]

            dest_dir = os.path.join(success_dir, acc)
            if not os.path.exists(dest_dir):
                os.mkdir(dest_dir)

#             shutil.copy(os.path.join(out_dir, out_pdb),
#                         os.path.join(dest_dir, out_pdb))

#             with open(os.path.join(dest_dir, acc + "_" + index + ".seq"), "wb") as seq_f:
#                 seq_f.write(">{acc}\n{seq}\n{struct}".format(acc=acc,
#                                                              seq=dict_digested_data[acc]["sequence"],
#                                                              struct=struct))

            with open(os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/for_pseudoviewer", acc + "_" + index + ".txt"), "wb") as seq_f:
                seq_f.write(">{acc}\n{seq}\n{struct}".format(acc=acc,
                                                             seq=dict_digested_data[acc]["sequence"],
                                                             struct=struct))

#             with open(os.path.join(dest_dir, acc + "_" + index + ".mat"), "wb") as mat_f:
#                 mat_f.write("\n".join([">{acc}\n{seq}".format(acc=dict_mat["accession"],
#                                                               seq=dict_mat["sequence"]) for \
#                                        dict_mat in dict_digested_data[acc]["matures"]]))