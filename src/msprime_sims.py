import msprime as ms
import tskit as ts

#need to rewrite as a function that just takes the species tree and spits out a newick

with open("/mnt/d/Projects/CoEvolution/data/dipterans/trees/astral_consensus.newick",'r') as file:
    data = file.read().replace('\n','')

dem = ms.Demography.from_species_tree(data,1000000,time_units="yr",generation_time=1000000000)
samples = {i.name:1 for i in dem.populations}

ts = ms.sim_ancestry(samples,demography=dem,sequence_length=1000)
tree = ts.first()
text_tree = tree.newick(precision=3)
for i in range(73,1):
    pop = dem.populations[i]
    name = pop.name
    text_tree=text_tree.replace(str(i)+":",name+":")

