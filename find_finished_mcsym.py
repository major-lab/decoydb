import os
import datetime

decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/decoys_100/"
completed_decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/completed_decoy_10_mar"

list_acc = os.listdir(decoy_dir)

for acc in sorted(list_acc):
    sub_decoy = os.listdir(os.path.join(decoy_dir, acc))

    now = datetime.datetime.now()
    list_finished = [sdecoy for sdecoy in sub_decoy if \
                     os.path.exists(os.path.join(decoy_dir, acc, sdecoy, "structure-report.txt")) or \
                     (datetime.datetime.fromtimestamp(os.path.getmtime(os.path.join(decoy_dir, acc, sdecoy))) <= now-datetime.timedelta(hours=2))]

#     if len(sub_decoy) != len(list_finished):
#         print len(sub_decoy), len(list_finished)

    if list_finished:
        if not os.path.exists(os.path.join(completed_decoy_dir, acc)):
            print "mkdir " + os.path.join(completed_decoy_dir, acc)
        print "\n".join(["mv {a} {b}".format(a=os.path.join(decoy_dir, acc, sdecoy),
                                             b=os.path.join(completed_decoy_dir, acc, sdecoy)) for sdecoy in list_finished])
