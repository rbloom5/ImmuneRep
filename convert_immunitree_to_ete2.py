import numpy
from ete2 import Tree
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#This function takes as an input the "node-based" tree repersentation output from ImmuniTree
#and outputs an ete2 Tree object

def convert_from_Immunitree(filelocation):

    DataIn = numpy.loadtxt(filelocation, dtype='string', delimiter=', ', skiprows=1 )

    node_dict = {}

    #cycles through each node in the .txt file and creates a leaf
    for row in DataIn:
    	#names the node
        name_string = 'node_'+ row[0]

        #Finds the root node and makes that the tree object
        if row[1] == '0':
            node_dict[name_string] = Tree(name=name_string)
        else:
        	#creates a string with the parent name
            parent = node_dict['node_'+ row[1]]
            #Creates a child with the node name under the node parent and saves it in a dict
            node_dict[name_string] = parent.add_child(name=name_string)

        #Adds the feature output from ImmuniTree to the node objects
        node_dict[name_string].add_features(node_size=int(row[2]), mutations=row[4], sequence=Seq(row[5], generic_dna))

    #returns the root node (which repersents the tree)
    return node_dict['node_1']
    #print tree.get_ascii(show_internal=True)