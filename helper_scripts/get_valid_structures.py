import os

decoy_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/processed/decoy/"

flashfold_dir = "/u/leongs/reproduction_projet_naim/rel20/2D/filtered"

list_acc = os.listdir(flashfold_dir)

for acc in sorted(list_acc):
    list_struct = []
    with open(os.path.join(flashfold_dir, acc), 'rb') as ff_o:
        list_struct = [elem.split()[0].strip() for elem in ff_o]

    valid_dir = os.path.join(decoy_dir, acc, 'out')
    print len(os.listdir(valid_dir))