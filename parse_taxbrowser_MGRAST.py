#!/usr/bin/env python

# this file uses the NCBI taxonomy to return the species, genus, family, order, class, phylum and kingdom of a taxon
# for multiple taxa a least common ancestor is returned
# 4 inputs: - the path to a folder containing the NCBI taxdump
#           - a file with multiple annotations for the same gene
#           - a file with single annotations for genes
#           - the name of the output
# 1 output - a table for all genes in the same order as in the two input files with LCA, species, genus, family, order, class, phylum, kingdom (tab separated)

# the code is mostly borrowed from ???, with the exception of the assembly of the output table; Anna Heintz-Buschart, April 2014

import os
import sys

taxdir = sys.argv[1]
multifile = sys.argv[2]
singlefile = sys.argv[3]
outfile = sys.argv[4]


# Definition of the class Node
class Node:
    """Node"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, number if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if name_object.has_key(tax_id):
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is the root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list

# Function to find common ancestor between two nodes or more
def common_ancestor(node_list):
    global name_object
    list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
    for node in node_list:
        list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
        ancestral_list = []                             
        for i in list1:
            if i in list2:                         # Identify common nodes between the two genealogy
                ancestral_list.append(i)                 
        list1 = ancestral_list                     # Reassing ancestral_list to list 1.
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestra_list is the common ancestor of all nodes.
    return common_ancestor                         # Return a node


#############################
#                           #
#   Read taxonomy files     #
#                           #
#############################

######################
# 
# Load names defintion

name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

# Load  NCBI names file
namepath = taxdir + "names.dmp"
name_file =  open(namepath,"r")
while 1:
    line = name_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split("|")
    tax_id, name = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_reverse[name] = tax_id
    if tab[3] == "scientific name":
             name_dict[tax_id] = name          # ... and load them into dictionary
name_file.close()

######################
# 
# Load kingdom definition

name_dict_kingdom = {}          # Initialise dictionary with TAX_ID:NAME

# Load kingdom file made in R
kingdompath = taxdir + "kingdom.txt"
kingdom_file =  open(kingdompath,"r")
while 1:
    line = kingdom_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, kingdom = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_kingdom[tax_id] = kingdom          # ... and load them into dictionary
kingdom_file.close()

######################
# 
# Load phylum definition

name_dict_phylum = {}          # Initialise dictionary with TAX_ID:NAME

# Load phylum file made in R
phylumpath = taxdir + "phylum.txt"
phylum_file =  open(phylumpath,"r")
while 1:
    line = phylum_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, phylum = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_phylum[tax_id] = phylum          # ... and load them into dictionary
phylum_file.close()

######################
# 
# Load class definition

name_dict_pclass = {}          # Initialise dictionary with TAX_ID:NAME

# Load class file made in R
pclasspath = taxdir + "class.txt"
pclass_file =  open(pclasspath,"r")
while 1:
    line = pclass_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, pclass = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_pclass[tax_id] = pclass          # ... and load them into dictionary
pclass_file.close()


######################
# 
# Load order definition

name_dict_order = {}          # Initialise dictionary with TAX_ID:NAME

# Load order file made in R
orderpath = taxdir + "order.txt"
order_file =  open(orderpath,"r")
while 1:
    line = order_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, order = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_order[tax_id] = order          # ... and load them into dictionary
order_file.close()

######################
# 
# Load family definition

name_dict_family = {}          # Initialise dictionary with TAX_ID:NAME

# Load family file made in R
familypath = taxdir + "family.txt"
family_file =  open(familypath,"r")
while 1:
    line = family_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, family = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_family[tax_id] = family         # ... and load them into dictionary
family_file.close()

######################
# 
# Load genus definition

name_dict_genus = {}          # Initialise dictionary with TAX_ID:NAME

# Load genus file made in R
genuspath = taxdir + "genus.txt"
genus_file =  open(genuspath,"r")
while 1:
    line = genus_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, genus = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_genus[tax_id] = genus         # ... and load them into dictionary
genus_file.close()

######################
# 
# Load species definition

name_dict_species = {}          # Initialise dictionary with TAX_ID:NAME

# Load species file made in R
speciespath = taxdir + "species.txt"
species_file =  open(speciespath,"r")
while 1:
    line = species_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split(";")
    tax_id, species = tab[0], tab[1]       # Assign tax_id and name ...
    name_dict_species[tax_id] = species         # ... and load them into dictionary
species_file.close()

######################
# 
# Load taxonomy

# Define taxonomy variable
global name_object
name_object = {}


# Load taxonomy NCBI file
taxonomypath = taxdir + "nodes.dmp"
taxonomy_file = open(taxonomypath,"r")
while 1:
    line = taxonomy_file.readline()
    if line == "":
        break
    #print line
    line = line.replace("\t","")
    tab = line.split("|")
    
    tax_id = str(tab[0])
    tax_id_parent = str(tab[1])
    division = str(tab[4])

    # Define name of the taxid
    name = "unknown"
    if tax_id in name_dict:
        name = name_dict[tax_id]
    
    if not name_object.has_key(tax_id):
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id   = tax_id        # Assign tax_id
    name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
    name_object[tax_id].name     = name          # Assign name
    
    if  tax_id_parent in name_object:
        children = name_object[tax_id].children  # If parent is is already in the object
        children.append(tax_id)                  # ...we found its children.
        name_object[tax_id].children = children  # ... so add them to the parent
taxonomy_file.close()


##################################################################
#                                                                #
#  Common ancestors or organism with some phylogeny into file    #
#                                                                #
##################################################################

commonAncestor = []
multi_file =  open(multifile,"r")
while 1:
    line = multi_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    names = line.split(";")
    ids = []
    for i in names:
        if i in name_dict_reverse:
            ids.append(str(name_dict_reverse[i]))
        else:
            ids.append(str(1))
    lca = common_ancestor(ids)
    commonAncestor.append(lca)
multi_file.close()

singles = []
single_file =  open(singlefile,"r")
while 1:
    line = single_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    names = line.split(";")
    for i in names:
        if i in name_dict_reverse:
            singles.append(str(name_dict_reverse[i]))
        else:
            singles.append(str(1))
single_file.close()


out_file = open(outfile, "w")
out_file.write("LCA" + "\t"+ "species" +"\t"+ "genus" + "\t"+ "family" +"\t"+ "order" + "\t"+  "class" + "\t"+  "phylum" + "\t"+  "kingdom"+ "\n")
for item in commonAncestor:
    kingdom_name = phylum_name = pclass_name = order_name = family_name = genus_name = species_name = "unknown"
    for anc in name_object[item].genealogy():
        if anc in name_dict_kingdom:
            kingdom_name = name_dict_kingdom[anc]
        if anc in name_dict_phylum:
            phylum_name = name_dict_phylum[anc]
        if anc in name_dict_pclass:
            pclass_name = name_dict_pclass[anc]
        if anc in name_dict_order:
            order_name = name_dict_order[anc]
        if anc in name_dict_family:
            family_name = name_dict_family[anc]
        if anc in name_dict_genus:
            genus_name = name_dict_genus[anc]
        if anc in name_dict_species:
            species_name = name_dict_species[anc]
    out_file.write(name_dict[item] + "\t"+ species_name +"\t"+ genus_name + "\t"+ family_name +"\t"+ order_name + "\t"+ pclass_name + "\t"+ phylum_name + "\t"+  kingdom_name + "\n")
for item in singles:
    kingdom_name = phylum_name = pclass_name = order_name = family_name = genus_name = species_name = "unknown"
    for anc in name_object[item].genealogy():
        if anc in name_dict_kingdom:
            kingdom_name = name_dict_kingdom[anc]
        if anc in name_dict_phylum:
            phylum_name = name_dict_phylum[anc]
        if anc in name_dict_pclass:
            pclass_name = name_dict_pclass[anc]
        if anc in name_dict_order:
            order_name = name_dict_order[anc]
        if anc in name_dict_family:
            family_name = name_dict_family[anc]
        if anc in name_dict_genus:
            genus_name = name_dict_genus[anc]
        if anc in name_dict_species:
            species_name = name_dict_species[anc]
    out_file.write(name_dict[item] + "\t"+ species_name +"\t"+ genus_name + "\t"+ family_name +"\t"+ order_name + "\t"+ pclass_name + "\t"+ phylum_name + "\t"+  kingdom_name + "\n")
out_file.close()
