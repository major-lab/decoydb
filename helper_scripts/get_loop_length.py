import cPickle
import os


def find_loop_length_list(info_dict):
    best_dir = os.path.join("/u/leongs/reproduction_projet_naim/rel20/3D/best_struct_bis",
                            info_dict["accession"])

    list_loop = []
    if os.path.exists(best_dir):
        curr_dict = dict()
        with open(os.path.join(best_dir, info_dict["accession"]+".pk"), 'rb') as tpk:
            curr_dict = cPickle.load(tpk)
        curr_struct = curr_dict["stat_struct"]
        curr_seq = curr_dict["seq"]

        end_5p = max(curr_dict["range_complementary_5p"])
        start_3p = min(curr_dict["range_complementary_3p"])

        i = end_5p
        while i <= start_3p:
            if curr_struct[i] == "(" and curr_struct[i+1] == ".":
                j = i+1
                all_good = True
                while j <= start_3p:
                    if curr_struct[j] == "(":
                        all_good = False
                        break
                    if curr_struct[j+1] == ")":
                        break
                    j += 1
                if all_good:
                    list_loop.append(j-i)
            i += 1
        return list_loop, True, curr_struct
    else:
        return list_loop, False, ""




info_list = []

with open("/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk", 'rb') as pk:
    info_list = cPickle.load(pk)

list_res = ["mature\tloop\tseq\t\len", ]
with open("/u/leongs/demandes_patron/list_mir.txt", 'rb') as list_mir:
    for line in list_mir:
        name = "hsa-" + line.strip()
        name = "-".join(name.split("-")[0:3])
        found = False
        for info in info_list:
            list_matures = [elem["name"].split("-") for elem in info["matures"]]
            for mat in list_matures:
                mat_name = "-".join(mat[0:3])
                if mat_name.startswith(name):
                    if mat_name != name:
                        diff = mat_name[len(name):]
                        is_digit = False
                        for c in diff:
                            if diff.isdigit():
                                is_digit = True
                                break
                        if not is_digit:
                            list_loop, curr_found, struct = find_loop_length_list(info)
                            if curr_found:
                                found = True
                                res = "\t".join((name, info["name"], info["sequence"], struct, ",".join([str(el) for el in list_loop])))
                                if not res in list_res:
                                    list_res.append(res)
                    elif mat_name == name:
                        list_loop, curr_found, struct = find_loop_length_list(info)
                        if curr_found:
                            found = True
                            res = "\t".join((name, info["name"], info["sequence"], struct, ",".join([str(el) for el in list_loop])))
                            if not res in list_res:
                                list_res.append(res)
        if not found:
            print line.strip()
# print "\n".join(list_res)
