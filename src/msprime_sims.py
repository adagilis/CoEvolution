import msprime as ms
import tskit as ts
import numpy as np
#need to rewrite as a function that just takes the species tree and spits out a newick

with open("/mnt/d/Projects/CoEvolution/data/dipterans/trees/msprime_rescaled.newick",'r') as file:
    data = file.read().replace('\n','')


Ne= 10000000 #Ok, so it insists the tree is not ultrametric. This gets around it by making the branch lengths _real_ long, but I hate it. Consider non-functional until fixed. Nominally - the tree is already in coalescent units. So, should rescale by *2Ne to get the generation time. Ne in flies is about 10^8, so this is about what we are doing, just using the yr-> gen conversion to imitate it.
dem = ms.Demography.from_species_tree(data,Ne,time_units="gen")
samples = {i.name:1 for i in dem.populations[1:73]}

ts = ms.sim_ancestry(samples,demography=dem,sequence_length=1000)
tree = ts.first()
ids = {x:dem.populations[x].name for x in range(1,73)}
text_tree = tree.as_newick(precision=3,node_labels=ids)

