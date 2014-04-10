from lxml import etree
import re
import os
import argparse
import cPickle
import math


def recursive_delete_circle(elem):
    list_to_recurse = []
    desired_node = None
    for child in elem.getchildren():
        if child.attrib.get("fill", "none") in ["gold", "PaleGoldenrod", "LightSteelBlue"]:
            desired_node = elem
            elem.remove(child)

        list_to_recurse.append(child)

    if desired_node is None:
        for child in list_to_recurse:
            desired_node = recursive_delete_circle(child)
            if desired_node is not None:
                break
    return desired_node

def recursive_delete_nucleotide_text_nodes(node):
    curr_list_text_att = []
    for child in node.getchildren():
        if "text" in child.tag and "onmouseover" in child.attrib:
            temp_dict = dict(value=child.text)
            temp_dict.update(child.attrib)
            curr_list_text_att.append(temp_dict)
            node.remove(child)
        else:
            curr_list_text_att.extend(recursive_delete_nucleotide_text_nodes(child))
    return curr_list_text_att


def convert_to_color_shade(val):
    level = math.ceil(val*255)
    return int(level)


def reformat_svg(svg_filepath, out_filepath, dict_color=dict(), pk_dict=dict(), list_color_lvl=[], add_statistics=True):
    new_svg_content = ""
    with open(svg_filepath, 'rb') as svg_c:
        new_svg_content = svg_c.read()#.replace('font-family="Tahoma"', 'font-family="courier new"')\
                                      #.replace('font-size="9.8"', 'font-size="9.0"')

    with open(out_filepath, 'wb') as svg_w:
        svg_w.write(new_svg_content)

    doc = etree.parse(out_filepath)

    root_elem = doc.getroot()
    desired_node = recursive_delete_circle(root_elem)
    list_text_att = recursive_delete_nucleotide_text_nodes(root_elem)

    for elem in list_text_att:
        splitted = elem["onmouseover"].split("'")
        pattern = re.search(r'\A[0-9]+', splitted[1])
        if pattern:
            st = pattern.group()
        else:
            alt_pattern = re.search(r'\([0-9]+\)', splitted[1].split()[1])
            if alt_pattern:
                st = alt_pattern.group().replace("(", "").replace(")", "")

        dict_color[str(list_color_lvl[int(st)-1])].append(dict(subnode=elem,
                                                               index=int(st)-1))
        if add_statistics:
            elem["onmouseover"] = elem["onmouseover"].replace("')", ", {perc}%')".format(perc=int(pk_dict["list_stats"][int(st)-1]["not_paired"]*100)))

    for color, list_subnodes in dict_color.iteritems():
        color_node = etree.SubElement(desired_node,
                                      "g",
                                      fill="rgb(255,{diff},{diff})".format(diff=255-int(color)))
        for subnode_dict in list_subnodes:
            subnode_index = subnode_dict["index"]
            subnode = subnode_dict["subnode"]
            
            if subnode_index in pk_dict["mature_range"]:
                fill_dict = {"stroke": "Blue",
                             "stroke-width": "1"}
            else:
                fill_dict = {"stroke": "none"}
            
            g_node = etree.SubElement(color_node,
                                      "g",
                                      **fill_dict)
            circle_node = etree.SubElement(g_node,
                                           "circle",
                                           cx=str(float(subnode["x"])+2.82),
                                           cy=str(float(subnode["y"])-3.92),
                                           r="4.8")

    temp_parent_node = etree.SubElement(desired_node,
                                        "g",
                                        fill="black",
                                        stroke="none")
    for elem in list_text_att:
        subelem = etree.SubElement(temp_parent_node,
                                   "text",
                                   onmouseover=elem["onmouseover"],
                                   onmouseout=elem["onmouseout"],
                                   x=elem["x"],
                                   y=elem["y"])
        subelem.text = elem["value"]

    outfile = open(out_filepath, 'w')
    doc.write(outfile)


def process_precursor(pk_filepath, svg_filepath, out_path, pv_dir, accession):
    pickled_dict = dict()
    with open(pk_filepath, 'rb') as pk_f:
        pickled_dict = cPickle.load(pk_f)

    list_color_lvl = [convert_to_color_shade(elem["not_paired"]) for elem in pickled_dict.get("list_stats")]
    dict_color = dict((str(elem), []) for elem in set(list_color_lvl))

    # first do the merged one
    reformat_svg(svg_filepath,
                 os.path.join(out_path, accession + ".svg"),
                 dict_color=dict_color,
                 pk_dict=pickled_dict,
                 list_color_lvl=list_color_lvl,
                 add_statistics=True)

    # now do the individual ones
    for individual in [elem for elem in os.listdir(pv_dir) if elem.startswith(accession) and \
                                                              "_" in elem and \
                                                              elem.endswith(".svg")]:
        list_color_lvl = [0 for elem in pickled_dict.get("list_stats")]
        dict_color = {"0":[]}
        reformat_svg(os.path.join(pv_dir, individual),
                     os.path.join(out_path, individual),
                     dict_color=dict_color,
                     pk_dict=pickled_dict,
                     list_color_lvl=list_color_lvl,
                     add_statistics=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--best_struct_dir', '-b', action="store", required=True, dest="best_struct_dir")
    parser.add_argument('--pseudoviewer_dir', '-p', action="store", required=True, dest="pseudoviewer_dir")
    parser.add_argument('--out_dir', '-o', action="store", required=True, dest="out_dir")

    ns = parser.parse_args()

    pseudoviewer_dir = ns.pseudoviewer_dir
    best_struct_dir = ns.best_struct_dir
    out_dir = ns.out_dir

    list_svg = sorted([elem for elem in os.listdir(pseudoviewer_dir) if elem.endswith(".svg") and not "_" in elem])

    for svg in sorted(list_svg):
        acc = svg.replace(".svg", "")
        pk_filepath = os.path.join(best_struct_dir, acc, acc + ".pk")
        svg_filepath = os.path.join(pseudoviewer_dir, svg)
        out_path = os.path.join(out_dir, acc)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        process_precursor(pk_filepath, svg_filepath, out_path, pseudoviewer_dir, acc)