#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "24/03/2014"
__version__ = "$Revision: 1.0"

from optparse import OptionParser
from operator import itemgetter
from collections import deque
#from igraph import *
import sys
import time
import os


#*****************************************************************************#
#                                                                             #
#                                Data structure                               #
#                                                                             #
#*****************************************************************************#
class Adjacency_list(object):
    neighbours = []
    parent = -1
    name = ''
    num = -1
    abundance = 0
    belongingRoot = -1
    checked = False
    seed = False

    # Initializer
    def __init__(self, neighbours):
        self.neighbours = neighbours


def make_adjlist(neighbours):
    adjacency_list = Adjacency_list(neighbours)
    return adjacency_list


#*****************************************************************************#
#                                                                             #
#                                   Helpers                                   #
#                                                                             #
#*****************************************************************************#
def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Break swarms."""

    parser = OptionParser(usage="usage: %prog --data_file filename",
                          description=desc,
                          version="%prog version 1.0")

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

    parser.add_option("-d", "--data_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="data_file",
                      help="set <FILENAME> as data file.")

    parser.add_option("-t", "--threshold",
                      metavar="<VALUE>",
                      type=float,
                      default=100,
                      action="store",
                      dest="threshold",
                      help="set <VALUE> as threshold.")

    parser.add_option("-m", "--manual",
                      default=False,
                      action="store_true",
                      dest="manualCut",
                      help="Activate manual cutting mode.")

    parser.add_option("-p", "--parameters",
                      default=True,
                      action="store_false",
                      dest="parameters",
                      help="Deactivate the use of parameters.")

    (options, args) = parser.parse_args()

    return options.fasta_file, options.swarm_file, options.data_file, \
        options.threshold, options.manualCut, options.parameters


def find_path(graph, start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    for node in graph[start].neighbours:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath:
                return newpath
    return None


### Two versions for rewire nodes' belonging root, if under treshold.
### Testing which is better - rewireFromRoot seems most promising.
def rewireFromNode(G, node, threshold):
    """
    Rewire node from breaking point, if belonging root is below threshold.
    i.e. performs BFS from breaking point to find new belonging root.
    """
    newParent = -1
    belongingRootAbundance = 0
    queue = deque([node])
    while newParent < 0:
        for neighbour in G[queue.popleft()].neighbours:
            curAbundance = G[G[neighbour].belongingRoot].abundance
            if curAbundance > threshold \
               and curAbundance > belongingRootAbundance:
                newParent = neighbour
                belongingRootAbundance = curAbundance
            else:
                queue.append(neighbour)
    #G[node].parent = newParent
    G[node].belongingRoot = G[newParent].belongingRoot


def rewireFromRoot(G, node, threshold):
    """
    Rewire node (at break point) from belonging root,
    if belonging root is below threshold.
    i.e. performs BFS from belonging root to find new belonging root.
    """
    newParent = -1
    belongingRootAbundance = 0
    queue = deque([G[node].belongingRoot])
    while newParent < 0:
        for neighbour in G[queue.popleft()].neighbours:
            curAbundance = G[G[neighbour].belongingRoot].abundance
            if curAbundance > threshold \
               and curAbundance > belongingRootAbundance:
                newParent = neighbour
                belongingRootAbundance = curAbundance
            else:
                queue.append(neighbour)
    #G[node].parent = newParent
    G[node].belongingRoot = G[newParent].belongingRoot


def outputSwarmFile(G, new_swarms, swarm_file):
    """
    Output new swarm file
    """
    tim = time.clock()

    output_file_swarm = os.path.splitext(swarm_file)[0]+"_new.swarm"
    with open(output_file_swarm, 'w') as f:
        for swarm in new_swarms:
            for node in swarm:
                f.write(str(G[node[0]].name)+"_"+str(node[1])+" ")
            f.write("\n")
    f.close()

    print "Time used to make output files:", time.clock()-tim


#*****************************************************************************#
#                                                                             #
#                                  Initializer                                #
#                                                                             #
#*****************************************************************************#
def buildGraph(fasta_file, swarm_file, data_file):
    """
    Set up data structure (graph)
    """
    print "Building data structure"
    tim = time.clock()

    with open(swarm_file, "rU") as swarm_file:
        for line in swarm_file:
            amplicons = [(amplicon.split("_")[0], int(amplicon.split("_")[1]))
                         for amplicon in line.strip().split(" ")]

    amplicon_index = {amplicon[0]: i for (i, amplicon) in enumerate(amplicons)}
    #nodeNames = [name[0] for name in amplicons]
    #abundance = [abund[1] for abund in amplicons]

    # Initialize graph structure
    G = [make_adjlist([]) for i in range(len(amplicon_index))]
    print "Network size:", len(G)

    # Insert name and abundance for each node
    for i in range(len(G)):
        G[i].name = amplicons[i][0]  # The nodes hashed name
        G[i].abundance = int(amplicons[i][1])  # the node's abundance
        G[i].num = i  # ID in graph

    if data_file:
        # Create list of neighbours
        with open(data_file, "rU") as data_file:
            for line in data_file:
                ampliconA, ampliconB, differences = line.split()
                ampliconA = amplicon_index[ampliconA]
                ampliconB = amplicon_index[ampliconB]
                G[ampliconA].neighbours.append(ampliconB)
                G[ampliconB].neighbours.append(ampliconA)
    #elif fasta_file:
        ### CODE FOR FASTA FILE HERE ###
    else:
        print "ERROR: NO DATA FILE OR FASTA FILE GIVEN"
        sys.exit(0)

    print "Time:", time.clock()-tim

    return G


#*****************************************************************************#
#                                                                             #
#                               Pretty printers                               #
#                                                                             #
#*****************************************************************************#
def manualCutter(G, possibleCuts):
    print
    print "ENTER MANUAL CUTTING MODE:"
    print "Instructions: To perform cut press y (yes) or enter, ",
    print "to keep edge press n (no)"
    print "Possible cuts:"
    print "[Cand for cut]    [Cand abundance]    ",
    print "[Closest root]    [Roots abundance]"
    space1 = 18
    space2 = 20
    space3 = 18
    finalCuts = []
    for edge in possibleCuts:
        # Ignore candidates if both belongs to same root.
        if G[edge[0]].belongingRoot != G[edge[1]].belongingRoot:
            len2 = space1-len(str([G[node].num for node in edge]))
            len3 = space2-len(str([G[node].abundance for node in edge]))-4
            len4 = space3-len(str([G[node].belongingRoot for node in edge]))-1
            print " ", [G[node].num for node in edge], " " * len2, \
                [G[node].abundance for node in edge], " " * len3, \
                [G[node].belongingRoot for node in edge], " " * len4, \
                [G[G[node].belongingRoot].abundance for node in edge],
            answer = raw_input("cut? ")
            if not answer == "n" or answer == "no" or answer == "nn":
                finalCuts.append(edge)
    return finalCuts


def prettyPrintCuts(G, finalCuts, tim):
    """
    A bad trial to give a nice output of the cuts.
    """
    print "[Cand for cut]    [Cand abundance]    ",
    print "[Closest root]    [Roots abundance]"
    space1 = 18
    space2 = 20
    space3 = 18
    for edge in finalCuts:
        len2 = space1-len(str([G[node].num for node in edge]))
        len3 = space2-len(str([G[node].abundance for node in edge]))-4
        len4 = space3-len(str([G[node].belongingRoot for node in edge]))-1
        print " ", [G[node].num for node in edge], " " * len2, \
            [G[node].abundance for node in edge], " " * len3, \
            [G[node].belongingRoot for node in edge], " " * len4, \
            [G[G[node].belongingRoot].abundance for node in edge]
        #print [G[node].name for node in edge]
    print "Number of final cuts:", len(finalCuts)
    print "Time:", tim, "\n"


#*****************************************************************************#
#                                                                             #
#                            Tools for break swarm                            #
#                                                                             #
#*****************************************************************************#
def assignParent(G, threshold):
    """
    Assign biggest neighbour as parent for each node in graph
    """
    print "\nAssigning parents"
    tim = time.clock()

    for node in G:
        # If node is a leaf, set only neighbour as parent if node is small
        if len(node.neighbours) == 1 and node.abundance < threshold:
            node.parent = node.neighbours[0]
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
                node.parent = -2  # Has no parent
                node.belongingRoot = node.num

    print "Time:", time.clock()-tim

    # Find name of node - Test purpose
    # if node.abundance == 2708:
    #     print "1:",node.num
    # if node.abundance == 3403:
    #     print "2:",node.num


def findPossibleCuts(G):
    """
    Find possible cuts by finding edges where no nodes is
    pointing from or to it.
    """
    print "\nFinding possible cuts...",
    tim = time.clock()

    possibleCuts = []
    for node in G:
        for neighbour in node.neighbours:
            if node.num < neighbour:
                if node.parent != neighbour and \
                   G[neighbour].parent != node.num:
                    possibleCuts.append((node.num, neighbour))
    #print "Possible cuts:\n",possibleCuts
    print "found", len(possibleCuts), "possible cuts"
    #prettyPrintCuts(G,possibleCuts,tim)

    print "Time:", time.clock() - tim
    return possibleCuts


def findBelongingRoot(G, possibleCuts):
    """
    Find corresponding root for each node in possible cuts
    """
    print "\nFinding belonging roots"
    tim = time.clock()

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

    print "Time:", time.clock()-tim


def rewireNode(G, possibleCuts, threshold):
    """
    "Rewire" nodes connected to a parent with belonging root
    that has abundance lower than threshold
    """
    print "\nRewiring nodes"
    tim = time.clock()

    for edge in possibleCuts:
        for node in edge:
            if G[G[node].belongingRoot].abundance < threshold:
                #rewireFromNode(G, node, threshold)
                rewireFromRoot(G, node, threshold)
    print "Time:", time.clock()-tim


def findFinalCuts(G, possibleCuts, threshold, manualCut, parameters):
    """
    Find final cuts, either by manually deciding the cuts or by
    using a parameter, or only using the threshold as tiebreaker.
    """
    print "\nFinding final cuts:"
    ### Measure time ###
    tim = time.clock()
    
    finalCuts = []
    # If manual cut on, the user gets to decide which edges will be cut.
    if manualCut:
        finalCuts = manualCutter(G, possibleCuts)
    # Automatic cut
    else:
        for edge in possibleCuts:
            # Ignore candidates if both belongs to same root.
            if G[edge[0]].belongingRoot != G[edge[1]].belongingRoot:
                ## TEST PURPOSE ONLY ##
                if parameters:
                    weakSpot = min(G[edge[0]].abundance,
                                   G[edge[1]].abundance)
                    biggestRoot = max(G[G[edge[0]].belongingRoot].abundance,
                                      G[G[edge[1]].belongingRoot].abundance)
                    smalletsRoot = min(G[G[edge[0]].belongingRoot].abundance,
                                       G[G[edge[1]].belongingRoot].abundance)
                    if (smalletsRoot/weakSpot > 25
                        and biggestRoot/smalletsRoot < 10) \
                            or smalletsRoot / weakSpot > threshold/2:
                        finalCuts.append(edge)
                else:
                    finalCuts.append(edge)

    ### Print Final Cuts ###
    if finalCuts:
        tim = time.clock() - tim
        prettyPrintCuts(G, finalCuts, tim)

    return finalCuts


def updateDataStructure(G, finalCuts):
    """
    Performing final cuts on graph structure and
    obtain possible seeds for new swarms.
    """
    tim = time.clock()
    new_swarm_seeds = [0]
    for edge in finalCuts:
        for node in edge:
            if not G[node].seed:
                new_swarm_seeds.append(node)
                G[node].seed = True
        G[edge[0]].neighbours.remove(edge[1])
        G[edge[1]].neighbours.remove(edge[0])
    print "Updating final cuts in data structure:", time.clock()-tim
    #print "Seeds:\n", new_swarm_seeds
    print "Number of possible seeds:", len(new_swarm_seeds), "\n"

    return new_swarm_seeds


#*****************************************************************************#
#                                                                             #
#                                 Break swarm                                 #
#                                                                             #
#*****************************************************************************#
def breakSwarm(G, threshold, manualCut, parameters):
    """
    Compute final cuts in the graph.
    """
    ### Set THRESHOLD value ###
    # Default: Ignore roots below 100.
    # if threshold >= 1, then we ignore roots below threshold,
    # else we ignore a percentage of the whole network size.
    threshold = threshold if threshold >= 1 else threshold*len(G)

    ### Assign a parent (biggest neighbour) to each node ###
    assignParent(G, threshold)

    ### Find possible cuts ###
    possibleCuts = findPossibleCuts(G)

    ### Find belonging root to each node in possibleCuts ###
    findBelongingRoot(G, possibleCuts)

    ### rewire nodes with small belonging roots ###
    rewireNode(G, possibleCuts, threshold)

    ### find final cuts: manually, parameter, or only by threshold ###
    finalCuts = findFinalCuts(G, possibleCuts, threshold, manualCut, parameters)

    # For testing - to see paths
    # print "Path:"
    # #path = find_path(G,12459,5889)
    # path = find_path(G, 0, 824)
    # print path
    # print [G[i].abundance for i in path]
    # print [G[i].belongingRoot for i in path]
    # print [G[G[i].belongingRoot].abundance for i in path]
    # #print [G[i].parent for i in path]
    # #print [G[i].name for i in path]
    # #print [G[i].neighbours for i in path]
    # print
    # nod = 1683
    # print G[nod].abundance
    # print G[nod].belongingRoot
    # print G[G[nod].belongingRoot].abundance
    # print G[nod].neighbours
    # print [G[G[nod].neighbours[i]].abundance
    #        for i in range(len(G[nod].neighbours))]

    ### Update data structure with final cuts
    new_swarm_seeds = updateDataStructure(G, finalCuts)

    return new_swarm_seeds


def findNewSwarms(G, seeds):
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
                new_swarm.append((currentNode, G[currentNode].abundance))
                for node in G[currentNode].neighbours:
                    if not G[node].checked:
                        queue.append(node)
                        G[node].checked = True

            new_swarm.sort(key=itemgetter(1), reverse=True)
            new_swarms.append(new_swarm)
    print "Visited nodes:", count
    print "BFS time:", time.clock()-tim
    print "Num swarms:", len(new_swarms)
    if count != len(G):
        print "ERROR: DATA LOST AFTER CUTTING."
    return new_swarms


#*****************************************************************************#
#                                                                             #
#                                     Main                                    #
#                                                                             #
#*****************************************************************************#
def main():

    totim = time.clock()

    ### Parse command line options ###
    fasta_file, swarm_file, data_file, threshold, manualCut, parameters = option_parse()

    ### Build data structure ###
    G = buildGraph(fasta_file, swarm_file, data_file)

    ### Compute cuts and break swarm ###
    new_swarm_seeds = breakSwarm(G, threshold, manualCut, parameters)

    ### Find new swarms ###
    new_swarms = findNewSwarms(G, new_swarm_seeds)

    ### Output new swarm file ###
    outputSwarmFile(G, new_swarms, swarm_file)

    print "Total time used:", time.clock() - totim


if __name__ == '__main__':

    main()

    sys.exit(0)


"""
Run example:
python breakOTUs.py -d file.data -s file.swarm
"""
