import os
from matplotlib import pyplot as plt

tprs, plfs = [], []
for item in os.listdir('.'):
    if os.path.isdir(item):
        for subitem in os.listdir(item):
            _path = os.path.join(item, subitem)
            if os.path.isdir(_path) and subitem == "collection_pmf":
                tprs.append(os.path.join(_path, "collection_pmf.tpr"))
                plfs.append(os.path.join(_path, "collection_pmf_{f}_F.xvg"\
                    .format(f = str(item))))

print(tprs)
print(plfs)

with open("tpr.dat", "w") as tpr:
    for item in tprs:
        tpr.write("%s\n" % (item))

with open("plf.dat", "w") as plf:
    for item in plfs:
        plf.write("%s\n" % (item))

os.system("gmx wham -it tpr.dat -if plf.dat -b 0")