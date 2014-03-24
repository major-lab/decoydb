import os

filtered_dir = "/u/leongs/reproduction_projet_naim/rel20/2D/filtered_alternative"

for elem in sorted(os.listdir(filtered_dir)):
    with open(os.path.join(filtered_dir, elem), 'rb') as st:
        if len(st.readlines()) < 100:
            print elem