import os

completed_list = ["completed_decoy_10_mar", ]

threeD_dir = "/u/leongs/reproduction_projet_naim/rel20/3D"

processed_list = ["processed_10_mar", ]

dict_success = dict()
for proc in processed_list:
    decoy_dir = os.path.join(threeD_dir, proc, "decoy")
    for acc in os.listdir(decoy_dir):
        out_dir = os.path.join(decoy_dir, acc, "out")
        for succ in os.listdir(out_dir):
            index = succ.split("_")[1].split("-")[0]
            if not acc in dict_success:
                dict_success[acc] = []
            dict_success[acc].append(int(index))

dict_fail = dict()
for fail in os.listdir("/u/leongs/reproduction_projet_naim/rel20/3D/failed_processed"):
    acc = fail.split("_")[0]
    index = fail.split("_")[1].split("-")[0]

    if not acc in dict_fail:
        dict_fail[acc] = []
    dict_fail[acc].append(int(index))

print len(dict_fail)

for c in completed_list:
    dir_path = os.path.join(threeD_dir, c)
    list_acc = os.listdir(dir_path)
    for acc in list_acc:
        acc_dir = os.path.join(dir_path, acc)
        for sub_acc in os.listdir(acc_dir):
            index = int(sub_acc.split("_")[1])
            if acc in dict_success:
                if not index in dict_fail.get(acc, []) and not index in dict_success.get(acc, []):
                    print "rm -rf " + os.path.join(acc_dir, sub_acc)
    print 'echo "finished {c}"'.format(c=c)