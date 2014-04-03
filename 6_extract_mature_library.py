import os
import gzip
import argparse
import cPickle

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--best_struct_dir', '-b', action="store", required=True, dest="best_struct_dir")
    parser.add_argument('--out_dir', '-o', action="store", required=True, dest="out_dir")

    ns = parser.parse_args()

    best_struct_dir = ns.best_struct_dir
    out_dir = ns.out_dir

    for acc in sorted(os.listdir(best_struct_dir)):
        # get the library_range for the microRNA
        pk_filepath = os.path.join(best_struct_dir, acc, acc + ".pk")
        info_dict = dict()
        with open(pk_filepath, "rb") as pk:
            info_dict = cPickle.load(pk)
        range_complementary_5p = info_dict["range_complementary_5p"]
        range_complementary_3p = info_dict["range_complementary_3p"]
        library_range = range_complementary_5p+range_complementary_3p

        pdb_file = [elem for elem in os.listdir(os.path.join(best_struct_dir, acc)) if elem.endswith(".pdb.gz")][0]

        # read the pdb file
        gz = gzip.open(os.path.join(best_struct_dir, acc, pdb_file))
        pdb_content = gz.readlines()
        gz.close()
        kept_content = []
        for line in pdb_content:
            stripped = line.strip()
            if stripped:
                # only keep content if the residue sequence number is inside the library_range
                current_ind = int(stripped[22:26])-1
                if current_ind in library_range:
                    kept_content.append(stripped)

        # write the library
        lib_dir = os.path.join(out_dir, acc)
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)

        # write the pdb file
        gz = gzip.open(os.path.join(lib_dir, acc + ".pdb.gz"), "wb")
        gz.write("\n".join(kept_content))
        gz.close()

        # write the script command for mcsymizer
        positionning = sorted((min(range_complementary_5p)+1, max(range_complementary_5p)+1,
                               min(range_complementary_3p)+1, max(range_complementary_3p)+1))

        # write the positionning into a file
        with open(os.path.join(lib_dir, "positions.txt"), 'wb') as pos_out:
            pos_out.write(",".join([str(elem) for elem in positionning]))