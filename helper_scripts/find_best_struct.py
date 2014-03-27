import os
import cPickle

def extract_pairing(structure):
    dict_open = dict()
    dict_close = dict()

    list_opener = []
    for i, char in enumerate(structure):
        if char == "(":
            list_opener.append(i)
        elif char == ")":
            # we add 1 to correct the fact that the 1st nt's index is 1 and not 0
            opener = list_opener.pop()
            closer = i
            dict_open[str(opener)] = closer
            dict_close[str(closer)] = opener

    return dict_open, dict_close


flashfold_dir = "/u/leongs/reproduction_projet_naim/rel20/2D/filtered"
digested_data_pk = "/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk"
processed_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/finished_processing_26_mar"
best_struct_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/best_struct"

list_digested_data = []
with open(digested_data_pk, 'rb') as dd:
    list_digested_data = cPickle.load(dd)

for index, hairpin_dict in enumerate(list_digested_data):
    hairpin_name = hairpin_dict['name']
    hairpin_acc = hairpin_dict['accession']
    hairpin_seq = hairpin_dict['sequence']

    if not os.path.exists(os.path.join(processed_dir, hairpin_acc)):
        continue

    mat_5p = ""
    mat_3p = ""
    for mature in hairpin_dict["matures"]:
        if mature["header"].endswith("5p"):
            mat_5p = mature["sequence"]
        else:
            mat_3p = mature["sequence"]

    mcfold_output = os.path.join(flashfold_dir, hairpin_acc)

    # get the structures that "worked"
    list_indexes = sorted([int(elem.split("-")[0].split("_")[-1]) for elem in os.listdir(os.path.join(processed_dir, hairpin_acc))])

    list_struct = []
    list_mcfold_l = []
    with open(mcfold_output, 'rb') as mcfold_o:
        list_mcfold_l = [elem.strip() for elem in mcfold_o.readlines()]
        list_struct = [list_mcfold_l[i] for i in list_indexes]

    fivep_start = -1
    fivep_end = -1
    if mat_5p:
        fivep_start = hairpin_seq.find(mat_5p)
        if fivep_start >= 0:
            fivep_end = fivep_start + len(mat_5p)

    threep_start = -1
    threep_end = -1
    if mat_3p:
        threep_start = hairpin_seq.find(mat_3p)
        if threep_start >= 0:
            threep_end = threep_start + len(mat_3p)

    struct_dict = dict()
    for struct in list_struct:
        fivep_struct = ""
        threep_struct = ""

        if fivep_start >= 0:
            fivep_struct = struct[fivep_start:fivep_end]

        if threep_start >= 0:
            threep_struct = struct[threep_start:threep_end]

        score = 0
        for struct2 in list_struct:
            if struct != struct2:
                curr_fivep_struct = ""
                curr_threep_struct = ""

                if fivep_start >= 0:
                    curr_fivep_struct = struct2[fivep_start:fivep_end]

                    for i in xrange(len(fivep_struct)):
                        if fivep_struct[i] == curr_fivep_struct[i]:
                            score += 1

                if threep_start >= 0:
                    curr_threep_struct = struct2[threep_start:threep_end]

                    for j in xrange(len(threep_struct)):
                        if threep_struct[j] == curr_threep_struct[j]:
                            score += 1
        struct_dict[struct] = score

    max_score = max(struct_dict.values())

    list_best_struct = []
    for struct, sc in struct_dict.iteritems():
        if sc == max_score:
            list_best_struct.append(struct)

    # sort the struct by energy
    list_best_struct.sort(key=lambda x: float(x.strip().split()[1]))

    best_struct = list_best_struct[0]

    # find the best_struct's index
    best_index = list_mcfold_l.index(best_struct)

    if not os.path.exists(os.path.join(best_struct_dir, hairpin_acc)):
        os.makedirs(os.path.join(best_struct_dir, hairpin_acc))
    output_path = os.path.join(best_struct_dir, hairpin_acc, "{acc}.struct".format(acc=hairpin_acc))
    with open(output_path, 'w') as out_file:
        out_file.write(">{hairpin_acc}\n{seq}\n{struct}".format(hairpin_acc=hairpin_acc,
                                                                seq=hairpin_seq,
                                                                struct=best_struct))

    # get all the structures that fit the best structure
    open_dict, close_dict = extract_pairing(best_struct.split()[0])
    range_5p = []
    range_3p = []
    if fivep_start:
        range_5p = range(fivep_start, fivep_end)
    if threep_start:
        range_3p = range(threep_start, threep_end)

    list_good_struct = []
    for struct in list_struct:
        this_open, this_close = extract_pairing(struct.split()[0])
        valid = True
        for op in range_5p:
            if not str(op) in this_open or not str(op) in open_dict:
                continue
            if this_open[str(op)] != open_dict[str(op)]:
                valid = False
                break
        for cl in range_3p:
            if not str(cl) in this_close or not str(cl) in close_dict:
                continue
            if this_close[str(cl)] != close_dict[str(cl)]:
                valid = False
                break
        if valid:
            list_good_struct.append(struct)

    list_stats = []

    for ind in xrange(len(best_struct.split()[0].strip())):
        none = len([elem for elem in list_good_struct if elem[ind] == "."])
        opened = len([elem for elem in list_good_struct if elem[ind] == "("])
        closed = len([elem for elem in list_good_struct if elem[ind] == ")"])

        list_stats.append(dict(not_paired=float(none)/float(len(list_good_struct)),
                               opened=float(opened)/float(len(list_good_struct)),
                               closed=float(closed)/float(len(list_good_struct))))

    # build the structure for pseudoviewer
    stat_struct = [c for c in best_struct.split()[0].strip()]
    for ind, s in enumerate(list_stats):
        if s["not_paired"] > 0.5:
            if stat_struct[ind] == "(":
                stat_struct[ind] = "."
                if str(ind) in open_dict:
                    stat_struct[open_dict[str(ind)]] = "."
            elif stat_struct[ind] == ")":
                stat_struct[ind] = "."
                if str(ind) in close_dict:
                    stat_struct[close_dict[str(ind)]] = "."

    to_be_pickled = dict(accession=hairpin_acc,
                         seq=hairpin_seq,
                         best_structure=best_struct,
                         list_best_structure=list_best_struct,
                         list_valid_structure=list_struct,
                         list_structure_for_bpstats=list_good_struct,
                         list_stats=list_stats,
                         stat_struct="".join(stat_struct),
                         mature_range=range_5p+range_3p)

    with open(os.path.join(best_struct_dir, hairpin_acc,"{acc}.pk".format(acc=hairpin_acc)),
              'wb') as pickle_file:
        cPickle.dump(to_be_pickled , pickle_file, -1)


    # make the file to be used for pseudoviewer
    with open(os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/for_pseudoviewer", "{acc}.txt".format(acc=hairpin_acc)),
              'wb') as out_file:
        out_file.write(">{hairpin_acc}\n{seq}\n{struct}".format(hairpin_acc=hairpin_acc,
                                                                seq=hairpin_seq,
                                                                struct="".join(stat_struct)))