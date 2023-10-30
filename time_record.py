import re
import os
import csv

def get_runtime(filename):
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r"(\d+\.\d+) seconds", line)
            if match:
                return float(match.group(1))
    return None

folder = "/home/tim/mixmodel/R_tree_comparison/MAFFT-AUTO-Untrimmed/result/"

with open('running_time.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Locus_name", "Model", "Run_time"])

    file_list = os.listdir(folder)
    for file in file_list:
        if not os.path.isfile(os.path.join(folder, file)):
            locus_name = file

            single_file = os.path.join(folder, file, f"Single_{locus_name}.iqtree")
            if os.path.exists(single_file):
                single_runtime = get_runtime(single_file)
                writer.writerow([locus_name, "Single", single_runtime])

            mix_file = os.path.join(folder, file, f"Mix_{locus_name}.iqtree")
            if os.path.exists(mix_file):
                mix_runtime = get_runtime(mix_file)
                writer.writerow([locus_name, "Mixture", mix_runtime])





