import subprocess
import os
import shutil
import matplotlib.pyplot as plt
from Bio import Phylo

# cwd = os.getcwd()


# if not os.path.exists('raxmlng_out'):
#     os.mkdir('raxmlng_out')
# elif os.listdir('raxmlng_out'):
#     for filename in os.listdir('raxmlng_out'):
#         os.remove(os.path.join(cwd,'raxmlng_out',filename))
# os.chdir(os.path.join(cwd,'raxmlng_out'))

# raxml_cmd = f"{cwd}/raxml -f x -m GTRGAMMA -p 12345 -x 12345 -# 100 -s {cwd}/trimmed_data.fasta -n out"
# raxmlng_cmd = f"{cwd}/raxml-ng --msa {cwd}/S1.raxml.reduced.phy --model GTR+G --prefix S1"
# subprocess.run(raxml_cmd,shell=True)
# if os.path.exists('trimmed_data.fasta.reduced'):
#     os.remove('trimmed_data.fasta.reduced')

# # ##
# tree = Phylo.read("RAxML_bipartitionsBranchLabels.out", "newick")
# for clade in tree.get_terminals():
#     if clade.name == 'matk_matk_contig_k41_1722':
#         clade.color == 'red'
# Phylo.draw(tree)
# for clade in tree.find_clades():
#     print(clade.name,clade.confidence)
Phylo.convert('RAxML_bipartitions.out','newick','test_tree','phyloxml')
tree = Phylo.read('test_tree','phyloxml')
leafs = tree.get_terminals()
 
clades = tree.find_clades()
clade_1 = next(clades)
print(dir(clade_1))
print(clade_1.confidence)
# for clade in clades:
#     if clade.name == 'MG225367.1_Cryptotaenia_canadensis':
#         clade.color = 'red'
#         print(clade.color)

# label_colors = {"MG225367.1_Cryptotaenia_canadensis": "red"}


# Phylo.draw(tree,show_confidence=True,label_colors=label_colors,do_show=False)
# plt.savefig("my_tree.png")






