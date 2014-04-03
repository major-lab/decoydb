import os
import json
import cPickle

dict_to_be_json_by_acc = dict()
dict_to_be_json_by_name = dict()
list_all_json = []
list_mature_json = []
list_precursor_json = []
dict_pk = dict()

with open("/u/leongs/reproduction_projet_naim/rel20/2D/digested_data.pk", 'rb') as dig_pk:
    dict_pk = cPickle.load(dig_pk)

for infodict in dict_pk:
    # do the precursor first
    accession = infodict["accession"]
    name = infodict["name"]
    sequence = infodict["sequence"]
    list_mature_acc = []
    list_mature_name = []

    if not accession in os.listdir("/u/leongs/reproduction_projet_naim/rel20/3D/best_struct_bis/"):
        continue

    # do the mature
    for matdict in infodict["matures"]:
        acc = matdict["accession"]
        nam = matdict["name"]
        seq = matdict["sequence"]
        if not acc in dict_to_be_json_by_acc:
            dict_to_be_json_by_acc[acc] = dict(accession=acc,
                                               name=nam,
                                               sequence=seq,
                                               list_precursor = [])
        dict_to_be_json_by_acc[acc]["list_precursor"].append(accession)
        if not nam in dict_to_be_json_by_name:
            dict_to_be_json_by_name[nam] = dict(accession=acc,
                                                name=nam,
                                                list_precursor = [])

        dict_to_be_json_by_name[nam]["list_precursor"].append(name)
        list_mature_acc.append(acc)
        list_mature_name.append(nam)

    dict_to_be_json_by_acc[accession] = dict(accession=accession,
                                             name=name,
                                             sequence=sequence,
                                             list_mature=list_mature_acc)
    dict_to_be_json_by_name[name] = dict(accession=accession,
                                         name=name,
                                         sequence=sequence,
                                         list_mature=list_mature_name)

for acc, infodict in dict_to_be_json_by_acc.iteritems():
    # if list_mature is in the dict, it means it's a precursor
    if 'list_mature' in infodict:
        list_obj = list_precursor_json
    else:
        list_obj = list_mature_json

    list_obj.append(dict(accession=infodict["accession"],
                         name=infodict["name"],
                         sequence=infodict["sequence"]))


with open("/u/leongs/reproduction_projet_naim/rel20/2D/dict_by_name.json", 'wb') as jsn:
    json.dump(dict_to_be_json_by_name, jsn)

with open("/u/leongs/reproduction_projet_naim/rel20/2D/dict_by_accession.json", 'wb') as jsa:
    json.dump(dict_to_be_json_by_acc, jsa)

with open("/u/leongs/reproduction_projet_naim/rel20/2D/list_mature.json", 'wb') as jsm:
    json.dump(list_mature_json, jsm)

with open("/u/leongs/reproduction_projet_naim/rel20/2D/list_precursor.json", 'wb') as jsp:
    json.dump(list_precursor_json, jsp)

with open("/u/leongs/reproduction_projet_naim/rel20/2D/list_all.json", 'wb') as jsall:
    json.dump(dict(aaData=list_mature_json+list_precursor_json), jsall)