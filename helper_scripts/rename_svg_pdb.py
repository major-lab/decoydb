import os
import shutil

# rename for_pseudoviewer
pv_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/for_pseudoviewer_bis"
new_pv_dir = "/u/leongs/reproduction_projet_naim/rel20/3D/new_pseudoviewer"
for f in sorted(os.listdir(pv_dir)):
    if not "_" in f:
        # directly copy it
        shutil.copy(os.path.join(pv_dir, f),
                    os.path.join(new_pv_dir, f))
    else:
        ext = f.split(".")[1]
        ind = int(f.split(".")[0].split("_")[1])
        acc = f.split("_")[0]
        if ext == "svg":
            shutil.copy(os.path.join(pv_dir, f),
                        os.path.join(new_pv_dir, "{acc}_{i}.{ext}".format(acc=acc, i=ind+1, ext=ext)))
        else:
            curr_text = ""
            with open(os.path.join(pv_dir, f), "rb") as fast_f:
                curr_text = fast_f.read()
            new_text = curr_text.replace("{acc}_{i}".format(acc=acc, i=ind), "{acc}_{i}".format(acc=acc, i=ind+1))
            with open(os.path.join(new_pv_dir, "{acc}_{i}.{ext}".format(acc=acc, i=ind+1, ext=ext)), "wb") as new_fast:
                new_fast.write(new_text)

