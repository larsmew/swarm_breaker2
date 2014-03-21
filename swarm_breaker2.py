__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "21/03/2014"
__version__ = "$Revision: 1.0"

from optparse import OptionParser
from operator import itemgetter
from collections import deque
from igraph import *
import tempfile
import subprocess
import sys
import time
import os

manuelCut = False #Default: False
THRESHOLD = 0
PARAMETER = True

class Adjacency_list(object):
    neighbours = []
    parent = -1
    name = ''
    num = -1
    abundance = 0
    belongingRoot = -1
    checked = False
    seed = False
    
    newParent = -1

    # Initializer 
    def __init__(self, neighbours):
        self.neighbours = neighbours

def make_adjlist(neighbours):
    adjacency_list = Adjacency_list(neighbours)
    return adjacency_list

def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Cut swarms."""
    
    parser = OptionParser(usage="usage: %prog --input_file filename",
        description = desc,
        version = "%prog version 1.0")
        
    parser.add_option("-f", "--fasta_file",
        metavar="<FILENAME>",
        action="store",
        dest="fasta_file",
        help="set <FILENAME> as fasta file.")
        
    parser.add_option("-s", "--swarm_file",
        metavar="<FILENAME>",
        action="store",
        dest="swarm_file",
        help="set <FILENAME> as swarm file.")
                      
    parser.add_option("-i", "--input_file",
        metavar = "<FILENAME>",
        action = "store",
        dest = "input_file",
        help = "set <FILENAME> as input file.")

    (options, args) = parser.parse_args()
    
    return options.fasta_file, options.swarm_file, options.input_file

def fasta_parse(fasta_file):
    """
    List amplicon ids, abundances and sequences, make a list and a dictionary
    """
    with open(fasta_file, "rU") as fasta_file:
        all_amplicons = dict()
        abundance = 0
        amplicon = 0
        for line in fasta_file:
            if line.startswith(">"):
                amplicon, abundance = line.strip(">\n").split("_")
            else:
                sequence = line.strip()
                all_amplicons[amplicon] = (int(abundance), str(sequence))
        return all_amplicons    

def run_modified_swarm(all_amplicons, swarm):
    """
    Write temporary fasta files, apply the modified swarm and collect
    the graph data
    """
    threshold = 1
    swarm_modified = ["/home/fred/Science/Projects/Swarms/swarm/swarm", "-d", str(threshold)]
    with open(os.devnull, "w") as devnull:
        with tempfile.SpooledTemporaryFile() as tmp_fasta_file:
            with tempfile.SpooledTemporaryFile() as tmp_swarm_results:
                for amplicon, abundance in swarm:
                    sequence = all_amplicons[amplicon][1]
                    #print(">", amplicon, "_", str(abundance), "\n", sequence, sep="", file=tmp_fasta_file)
                tmp_fasta_file.seek(0)  # rewind to the begining of the file
                proc = subprocess.Popen(swarm_modified,
                                        stderr=tmp_swarm_results,
                                        stdout=devnull,
                                        stdin=tmp_fasta_file,
                                        close_fds=True)
                proc.wait()  # usefull or not?
                tmp_swarm_results.seek(0)  # rewind to the begining of the file
                graph_data = [line.strip().split("\t")[1:4]
                              for line in tmp_swarm_results
                              if line.startswith("@")]
                return graph_data

def manualCutter(G,possibleCuts):
    print
    print "ENTER MANUAL CUTTING MODE:"
    print "Instructions: To perform cut press y (yes) or enter, to keep edge press n (no)"
    print "Possible cuts:"
    print "[Cand for cut]    [Cand abundance]    [Closest root]    [Roots abundance]"
    space1 = 18
    space2 = 20
    space3 = 18
    finalCuts = []
    for edge in possibleCuts:
        len2 = space1-len(str([G[node].num for node in edge]))
        len3 = space2-len(str([G[node].abundance for node in edge]))-4
        len4 = space3-len(str([G[node].belongingRoot for node in edge]))-1
        print " ",[G[node].num for node in edge]," "*len2,[G[node].abundance for node in edge]," "*len3,[G[node].belongingRoot for node in edge]," "*len4,[G[G[node].belongingRoot].abundance for node in edge],
        answer = raw_input("cut? ")
        if not answer == "n" or answer == "no" or answer == "nn":
            finalCuts.append(edge)
    return finalCuts

def prettyPrintCuts(G,finalCuts,tim):
    '''
    A bad trail to give a nice output of the cuts.
    '''
    print
    print "Final cuts:"
    print "[Cand for cut]    [Cand abundance]    [Closest root]    [Roots abundance]"
    space1 = 18
    space2 = 20
    space3 = 18
    for edge in finalCuts:
        len2 = space1-len(str([G[node].num for node in edge]))
        len3 = space2-len(str([G[node].abundance for node in edge]))-4
        len4 = space3-len(str([G[node].belongingRoot for node in edge]))-1
        print " ",[G[node].num for node in edge]," "*len2,[G[node].abundance for node in edge]," "*len3,[G[node].belongingRoot for node in edge]," "*len4,[G[G[node].belongingRoot].abundance for node in edge]
        #print [G[node].name for node in edge]
    print "Number of final cuts:",len(finalCuts)
    print "Time:",tim
    print
    

# Rewire from breaking point
def rewireFromNode(G,node,threshold):
    newParent = -1
    belongingRootAbundance = 0
    queue = deque([node])
    while newParent < 0:
        for neighbour in G[queue.popleft()].neighbours:
            curAbundance = G[G[neighbour].belongingRoot].abundance
            if curAbundance > threshold and curAbundance > belongingRootAbundance:
                newParent = neighbour
                belongingRootAbundance = curAbundance
            else:
                queue.append(neighbour)
    G[node].parent = newParent
    G[node].belongingRoot = G[newParent].belongingRoot

def rewireFromRoot(G,node,threshold):
    newParent = -1
    belongingRootAbundance = 0
    queue = deque([G[node].belongingRoot])
    while newParent < 0:
        for neighbour in G[queue.popleft()].neighbours:
            curAbundance = G[G[neighbour].belongingRoot].abundance
            if curAbundance > threshold and curAbundance > belongingRootAbundance:
                newParent = neighbour
                belongingRootAbundance = curAbundance
            else:
                queue.append(neighbour)
    G[node].parent = newParent
    G[node].belongingRoot = G[newParent].belongingRoot


def rewire2(G,node,cutNode):
    print "rewire2"
    newParent = -1
    belongingRootAbundance = G[G[node].belongingRoot].abundance
    queue = deque([node])
    while queue and newParent < 0:
        for neighbour in G[queue.popleft()].neighbours:
            if neighbour != G[node].parent and neighbour != cutNode:
                curAbundance = G[G[neighbour].belongingRoot].abundance
                if curAbundance > belongingRootAbundance:
                    newParent = neighbour
                    belongingRootAbundance = curAbundance
                elif G[node].belongingRoot != G[neighbour].belongingRoot:
                    queue.append(neighbour)
    if newParent == -1:
        return False
    G[node].parent = newParent
    G[node].belongingRoot = G[newParent].belongingRoot
    return True

def simpleCut(G):
    
    # Assign a "parent" to each node
    for node in G:
        if node.abundance == 2708:
            print "1:",node.num
        if node.abundance == 3403:
            print "2:",node.num
        # If node is a leaf, set only neighbour as parent
        if len(node.neighbours) == 1  and node.abundance < THRESHOLD:
            node.parent = node.neighbours[0]
            #node.belongingRoot = node.neighbours[0]
        # If node is the biggest in graph, ignore it
        elif node.num == 0:
            node.parent = -2
            node.belongingRoot = node.num
            node.seed = True
        # Choose neighbour with highest abundance as parent
        else:
            biggestNeighbour = -1
            biggestNeighbourAbundance = 0
            for neighbour in node.neighbours:
                if G[neighbour].abundance > biggestNeighbourAbundance:
                    biggestNeighbour = neighbour
                    biggestNeighbourAbundance = G[neighbour].abundance
            if biggestNeighbourAbundance >= node.abundance:
                node.parent = biggestNeighbour
            else:
                node.parent = -2 # Has no parent
                node.belongingRoot = node.num
    
    tim = time.clock()
    
    # Find possible cuts (in O(m))
    possibleCuts = []
    for node in G:
        for neighbour in node.neighbours:
            if node.num < neighbour:
                if node.parent != neighbour and G[neighbour].parent != node.num:
                    possibleCuts.append((node.num,neighbour))
    #print "Possible cuts:\n",possibleCuts
    print "Number of possible cuts:",len(possibleCuts)
    #prettyPrintCuts(G,possibleCuts,tim)
    print "Time:",time.clock()-tim
    
     
    tim = time.clock()
    
    # Find corresponding root for each node in possible cuts
    for edge in possibleCuts:
        for node in edge:
            if G[node].belongingRoot == -1:
                currentNode = G[node]
                #path = [node]
                while G[currentNode.parent].belongingRoot == -1:
                    #path.append(currentNode.parent)
                    
                    # Skip pair of nodes pointing at each other
                    if G[currentNode.parent].parent == currentNode.num:
                        G[currentNode.parent].belongingRoot = currentNode.num
                        break
                    currentNode = G[currentNode.parent]
                #for n in path:
                #    G[n].belongingRoot = G[currentNode.parent].belongingRoot
                G[node].belongingRoot = G[currentNode.parent].belongingRoot

    
    # "Rewire" nodes connected to a parent with abundance lower than threshold
    threshold = THRESHOLD if THRESHOLD > 1 else THRESHOLD*len(G)
    for edge in possibleCuts:
        for node in edge:
            if G[G[node].belongingRoot].abundance < threshold:
                #rewireFromNode(G,node,threshold)
                rewireFromRoot(G,node,threshold)
    
    finalCuts = []
    # If manual cut on, the user gets to decide which edges will be cut.
    if manuelCut:
        finalCuts = manualCutter(G,possibleCuts,finalCuts)
    # Automatic cut
    else: 
        for edge in possibleCuts:
            if G[edge[0]].belongingRoot != G[edge[1]].belongingRoot:
                ## TEST PURPOSE ONLY ##
                if PARAMETER:
                    weakSpot = min(G[edge[0]].abundance, G[edge[1]].abundance)
                    biggestRoot = max(G[G[edge[0]].belongingRoot].abundance, G[G[edge[1]].belongingRoot].abundance)
                    smalletsRoot = min(G[G[edge[0]].belongingRoot].abundance, G[G[edge[1]].belongingRoot].abundance)
                    if (smalletsRoot/weakSpot > 25 and biggestRoot/smalletsRoot < 10) or smalletsRoot / weakSpot > THRESHOLD/2:
                        finalCuts.append(edge)
                else:
                    finalCuts.append(edge)
    
    
    ### Print Final Cuts ###
    if finalCuts:
        tim = time.clock()-tim
        prettyPrintCuts(G,finalCuts,tim)

    return finalCuts

def findNewSwarms(G,seeds):
    ### Performing BFS to find new swarms ###
    print "Performing BFS to discover new swarms..."
    tim = time.clock()
    new_swarms = []
    count = 0
    for seed in seeds:
        if not G[seed].checked:
            new_swarm = []
            G[seed].checked = True
            # BFS
            queue = deque([seed])
            while queue:
                count += 1
                currentNode = queue.popleft()
                new_swarm.append((currentNode,G[currentNode].abundance))
                for node in G[currentNode].neighbours:
                    if not G[node].checked:
                        queue.append(node)
                        G[node].checked = True
            
            new_swarm.sort(key=itemgetter(1), reverse=True)
            new_swarms.append(new_swarm)
    print "Visited nodes:",count
    print "BFS time:",time.clock()-tim
    print "Num swarms:",len(new_swarms)
    if count != len(G):
        print "ERROR: DATA LOST AFTER CUTTING."
    return new_swarms

def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    for node in graph[start].neighbours:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None
    
def outputGML(G,possibleCuts,output):
    with open(output, 'w') as f:
        # Setup
        f.write("graph\n[\n")
        f.write("\thierarchic\t1\n")
        f.write("\tlabel\t\"\"\n")
        f.write("\tdirected\t1\n")
        
        for node in G:
            f.write("\tnode\n")
            f.write("\t[\n")
            f.write("\t\tid\t"+str(node.num)+"\n")
            f.write("\t\tlabel\t\""+str(node.abundance)+"\"\n")
            f.write("\t]\n")
        
        for node in G:
            if node.parent != -2:
                f.write("\tedge\n")
                f.write("\t[\n")
                f.write("\t\tsource\t"+str(node.num)+"\n")
                f.write("\t\ttarget\t"+str(node.parent)+"\n")
                f.write("\t\tgraphics\n")
                f.write("\t\t[\n")
                f.write("\t\t\tfill\t\"#000000\"\n")
                f.write("\t\t\ttargetArrow\t\"delta\"\n")
                f.write("\t\t]\n")
                f.write("\t]\n")
        
        for edge in possibleCuts:
            f.write("\tedge\n")
            f.write("\t[\n")
            f.write("\t\tsource\t"+str(edge[0])+"\n")
            f.write("\t\ttarget\t"+str(edge[1])+"\n")
            f.write("\t\tgraphics\n")
            f.write("\t\t[\n")
            f.write("\t\t\tfill\t\"#FF0000\"\n")
            f.write("\t\t\tsourceArrow\t\"delta\"\n")
            f.write("\t\t\ttargetArrow\t\"delta\"\n")
            f.write("\t\t]\n")
            f.write("\t]\n")
        
        # Finish
        f.write("]")
    f.close

if __name__ == '__main__':
    
    # Set THRESHOLD value
    THRESHOLD = 100
    
    # Parse command line options.
    fasta_file, swarm_file, input_file = option_parse()
    print fasta_file, swarm_file,input_file
    
    output_file_gml = "test.gml"
    output_file_swarm = os.path.splitext(swarm_file)[0]+"_new.swarm"
    output_file_output = os.path.splitext(swarm_file)[0]+"_new.data"
    
    with open(swarm_file, "rU") as swarm_file:
        for line in swarm_file:
            amplicons = [(amplicon.split("_")[0], int(amplicon.split("_")[1]))
                     for amplicon in line.strip().split(" ")]
        
    amplicon_index = {amplicon[0]: i for (i, amplicon) in enumerate(amplicons)}
    #nodeNames = [name[0] for name in amplicons]
    #abundance = [abund[1] for abund in amplicons]
    
    # Initialize graph structure
    G = [make_adjlist([]) for i in range(len(amplicon_index))]
    print "Network size:",len(G)
    
    # Insert name and abundance for each node
    for i in range(len(G)):
        G[i].name = amplicons[i][0] # The nodes hashed name
        G[i].abundance = int(amplicons[i][1]) # the node's abundance
        G[i].num = i # ID in graph
    
    if input_file is not None:
        # Create list of neighbours
        with open(input_file, "rU") as input_file:
            for line in input_file:
                ampliconA, ampliconB, differences = line.split()
                ampliconA = amplicon_index[ampliconA]
                ampliconB = amplicon_index[ampliconB]
                G[ampliconA].neighbours.append(ampliconB)
                G[ampliconB].neighbours.append(ampliconA)
    elif fasta_file is not None:
        all_amplicons = fasta_parse(fasta_file)
    else:
        print "ERROR: NO INPUT FILE OR FASTA FILE GIVEN"
        sys.exit(0)
    
    ### Compute simple cut ###
    moreCuts = True
    new_swarm_seeds = [0]
    iteration = 0
    while moreCuts:
        iteration += 1
        print "\nRunning iteration "+str(iteration)+":"
        finalCuts = simpleCut(G)
        if len(finalCuts) == 0:
            moreCuts = False
            print "No more final cuts found\n"
        else:

            ### Output GML file with possibe cut in red ###
            # tim = time.clock()
            # outputGML(G,finalCuts,output_file_gml)
            # print "Time to output GML:",time.clock()-tim
        
            if iteration == 1:
                print "Path:"
                #path = find_path(G,12459,5889)
                path = find_path(G,0,824)
                print path
                print [G[i].abundance for i in path]
                print [G[i].belongingRoot for i in path]
                print [G[G[i].belongingRoot].abundance for i in path]
                #print [G[i].parent for i in path]
                #print [G[i].name for i in path]
                #print [G[i].neighbours for i in path]
                print
                nod = 1683
                print G[nod].abundance
                print G[nod].belongingRoot
                print G[G[nod].belongingRoot].abundance
                print G[nod].neighbours
                print [G[G[nod].neighbours[i]].abundance for i in range(len(G[nod].neighbours))]
            

            ### Performing final cuts on graph structure ###
            tim = time.clock()
            for edge in finalCuts:
                for node in edge:
                    if not G[node].seed:
                        new_swarm_seeds.append(node)
                        G[node].seed = True
                G[edge[0]].neighbours.remove(edge[1])
                G[edge[1]].neighbours.remove(edge[0])
            print "Making final cuts:",time.clock()-tim
    
    # path = find_path(G,0,1)
    # print path
    # print [G[i].abundance for i in path]
    # print G[1683].neighbours
    # print [G[G[1683].neighbours[i]].abundance for i in range(len(G[1683].neighbours))]
    # print G[G[1683].belongingRoot].abundance
    # print G[1683].belongingRoot
    # print
    
    print "Seeds:\n",new_swarm_seeds
    print "Num possible seeds:",len(new_swarm_seeds)
    print 

    ### Find new swarms ###
    new_swarms = findNewSwarms(G,new_swarm_seeds)

    # print
    # print G[29689].belongingRoot
    # for seed in new_swarm_seeds:
    #     path = find_path(G,29689,seed)
    #     if path:
    #         print path
    
    tim = time.clock()
    ### Output new swarm file ###
    with open(output_file_swarm, 'w') as f:
        for swarm in new_swarms:
            for node in swarm:
                f.write(str(G[node[0]].name)+"_"+str(node[1])+" ")
            f.write("\n")
    f.close()
    
    ### Output new data file ###
    # with open(output_file_output, 'w') as f:
    #     for swarm in new_swarms:
    #         for node in swarm:
    #             for neighbour in G[node[0]].neighbours:
    #                 if node[1] > G[neighbour].abundance:
    #                     f.write(str(G[node[0]].name)+"\t"+str(G[neighbour].name)+"\t1\n")
    #                 elif node[1] == G[neighbour].abundance:
    #                     f.write(str(G[node[0]].name)+"\t"+str(G[neighbour].name)+"\t1\n")
    #                     G[neighbour].neighbours.remove(node[0])
    # f.close()
    print "Time to make output files:",time.clock()-tim
    
    # Pause program
    #raw_input()
    
    # i = 1
    # for swarm in new_swarms:
    #     if len(swarm) > 1500:
    #         outputGML(swarm,[],"test"+str(i)+".gml")
    #         i += 1       
            
    # print 
    # print "573:",G[573].parent
    # print "573:",G[573].abundance
    # print G[573].neighbours
    # print [G[neighbour].abundance for neighbour in G[573].neighbours]
    # print G[573].belongingRoot
    
    
    # paths = find_all_paths(G,0,1)
    # print paths
    # for path in paths:
    #     print [G[i].abundance for i in path]
    
    


"""
python -m cProfile breakOTUs.py -i OTU_001_26d2046e1a0f79450d6233ef29aab44e.data -s OTU_001_26d2046e1a0f79450d6233ef29aab44e.swarm
"""
