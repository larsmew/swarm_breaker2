#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

__author__ = "Lars Andersen <larsmew@gmail.com>"
__date__ = "30/04/2014"
__version__ = "$Revision: 1.2"

from optparse import OptionParser
from operator import itemgetter
from collections import deque
# from igraph import *
import tempfile
import subprocess
import sys
import time
import os


# *************************************************************************** #
#                                                                             #
#                                Data structure                               #
#                                                                             #
# *************************************************************************** #

class Node(object):
    """
    Define the data structure of a node and its direct neighbours
    """
    neighbours = []
    parent = -1
    name = ''
    num = -1
    abundance = 0
    belongingRoot = -1
    checked = False
    seed = False

    # Initializer
    def __init__(self):
        self.neighbours = []


def create_node():
    """
    Create node objects
    """
    node = Node()
    return node


# *************************************************************************** #
#                                                                             #
#                                   Helpers                                   #
#                                                                             #
# *************************************************************************** #

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

    parser.add_option("-b", "--binary_swarm",
                      metavar="<FILENAME>",
                      action="store",
                      dest="swarm_path",
                      help="path to binary swarm file.")

    parser.add_option("-m", "--manual",
                      default=False,
                      action="store_true",
                      dest="manualCut",
                      help="Activate manual cutting mode.")

    parser.add_option("-p", "--parameters",
                      default=False,
                      action="store_true",
                      dest="parameters",
                      help="Activate the use of parameters.")

    (options, args) = parser.parse_args()

    return options.fasta_file, options.swarm_file, options.data_file, \
        options.threshold, options.swarm_path, options.manualCut, \
        options.parameters


def find_path(graph, start, end, path=[]):
    """
    Find the shortest path between two nodes
    """
    path = path + [start]
    if start == end:
        return path
    for node in graph[start].neighbours:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath:
                return newpath
    return None


def outputSwarmFile(G, new_swarms, swarm_file, fasta_file):
    """
    Output new swarm file
    """
    tim = time.clock()

    # if swarm_file:
    #     output_file_swarm = os.path.splitext(swarm_file)[0] + "_new.swarms"
    # else:
    #     output_file_swarm = os.path.splitext(fasta_file)[0] + "_new.swarms"
    # with open(output_file_swarm, 'w') as f:
    for swarm in new_swarms:
        # print(" ".join([str(G[node[0]].name)+"_"+str(node[1])
        #      for node in swarm]), file=f)
        print(" ".join([node[2] + "_" + str(node[1])
              for node in swarm]), file=sys.stdout)

    print("Time used to make output files:", time.clock() - tim,
          file=sys.stderr)
    return None


# *************************************************************************** #
#                                                                             #
#                                  Initializer                                #
#                                                                             #
# *************************************************************************** #

def fasta_parse(fasta_file):
    """
    List amplicon ids, abundances and sequences, make a list and a dictionary
    """
    with open(fasta_file, "rU") as fasta_file:
        all_amplicons = dict()
        for line in fasta_file:
            if line.startswith(">"):
                amplicon, abundance = line.strip(">\n").rsplit("_", 1)
            else:
                sequence = line.strip()
                all_amplicons[amplicon] = (int(abundance), sequence)
        return all_amplicons


def swarm_parse(swarm_file):
    """
    List amplicons contained in each swarms, sort by decreasing
    abundance. Sort the list of swarms by decreasing mass and
    decreasing size.
    """
    with open(swarm_file, "rU") as swarm_file:
        swarms = list()
        for line in swarm_file:
            amplicons = [(amplicon.rsplit("_", 1)[0],
                         int(amplicon.rsplit("_", 1)[1]))
                         for amplicon in line.strip().split(" ")]
            # Sort amplicons by decreasing abundance and alphabetical order
            amplicons.sort(key=itemgetter(1, 0), reverse=True)
            top_amplicon, top_abundance = amplicons[0]
            swarm_size = len(amplicons)
            swarm_mass = sum([amplicon[1] for amplicon in amplicons])
            swarms.append([top_amplicon, swarm_mass, swarm_size,
                           top_abundance, amplicons])
        # Sort swarms on mass, size and seed name
        swarms.sort(key=itemgetter(1, 2, 0), reverse=True)
        return swarms


def run_swarm(binary, all_amplicons, swarm, threshold):
    """
    Write temporary fasta files, run swarm and collect the graph data
    """
    # swarm_command = [binary, "-b", "-d", str(threshold)]
    swarm_command = [binary, "-b", "-d", "1"]
    with open(os.devnull, "w") as devnull:
        with tempfile.SpooledTemporaryFile() as tmp_fasta_file:
            with tempfile.SpooledTemporaryFile() as tmp_swarm_results:
                for amplicon, abundance in swarm:
                    sequence = all_amplicons[amplicon][1]
                    print(">", amplicon, "_", str(abundance), "\n", sequence,
                          sep="", file=tmp_fasta_file)
                tmp_fasta_file.seek(0)  # rewind to the begining of the file
                try:
                    proc = subprocess.Popen(swarm_command,
                                            stderr=tmp_swarm_results,
                                            stdout=devnull,
                                            stdin=tmp_fasta_file,
                                            close_fds=True)
                except (OSError, 2):
                    print("Error ",
                          "Cannot find the swarm binary in the $PATH" +
                          "directories\n",
                          "Use -b /path/to/swarm to indicate the correct path",
                          sep="", file=sys.stderr)
                    sys.exit(-1)
                proc.wait()  # usefull or not?
                tmp_swarm_results.seek(0)  # rewind to the begining of the file
                graph_data = [line.strip().split("\t")[1:4]
                              for line in tmp_swarm_results
                              # if line.startswith("@")]
                              if "@@" in line]
                return graph_data


def build_graph(graph_data, amplicons, is_data_file):
    """
    List pairwise relations in a swarm. Note that not all pairwise
    relations are stored. That's why the graph exploration must always
    start from the most abundant amplicon, and must be reiterated for
    sub-swarms after a breaking.
    """
    amplicon_index = {amplicon[0]: i for (i, amplicon)
                      in enumerate(amplicons)}

    G = [create_node() for _ in range(len(amplicons))]

    for i in range(len(amplicon_index)):
        G[i].name = amplicons[i][0]  # The nodes hashed name
        G[i].abundance = int(amplicons[i][1])  # the node's abundance
        G[i].num = i  # ID in graph

    if is_data_file:
        with open(graph_data, "rU") as graph_data:
            for line in graph_data:
                ampliconA, ampliconB, differences = line.split()
                ampliconA = amplicon_index[ampliconA]
                ampliconB = amplicon_index[ampliconB]
                G[ampliconA].neighbours.append(ampliconB)
                G[ampliconB].neighbours.append(ampliconA)
    else:
        for line in graph_data:
            ampliconA, ampliconB, differences = line
            ampliconA = amplicon_index[ampliconA]
            ampliconB = amplicon_index[ampliconB]
            G[ampliconA].neighbours.append(ampliconB)
            G[ampliconB].neighbours.append(ampliconA)

    print("Graph size:", len(G), file=sys.stderr)

    return G


def initGraph(swarm_file, graph_data, is_data_file):
    """
    Initialize the graph structure and insert data
    """
    tim2 = time.clock()
    # Read amplicons from swarm file
    amplicons = []
    for line in swarm_file:
        amplicons += [(amplicon.rsplit("_", 1)[0],
                      int(amplicon.rsplit("_", 1)[1])) for
                      amplicon in line.strip().split(" ")]

    amplicon_index = {amplicon[0]: i for (i, amplicon)
                      in enumerate(amplicons)}
    print("Time for reading swarm file:", time.clock() - tim2)

    tim2 = time.clock()
    # Initialize graph structure
    G = [create_node() for _ in range(len(amplicon_index))]
    print("Time for creating graph:", time.clock() - tim2)

    tim2 = time.clock()
    # Insert name and abundance of each node in the graph
    for i in range(len(amplicon_index)):
        G[i].name = amplicons[i][0]  # The nodes hashed name
        G[i].abundance = int(amplicons[i][1])  # the node's abundance
        G[i].num = i  # ID in graph
    print("Time for inserting name, abundance, id:", time.clock() - tim2)

    tim2 = time.clock()
    for line in graph_data:
        if is_data_file:
            ampliconA, ampliconB, differences = line.split()
        else:
            ampliconA, ampliconB = line[0], line[1]
        ampliconA = amplicon_index[ampliconA]
        ampliconB = amplicon_index[ampliconB]
        G[ampliconA].neighbours.append(ampliconB)
        G[ampliconB].neighbours.append(ampliconA)
    print("Time for inserting nodes:", time.clock() - tim2)

    return G


def buildGraph(fasta_file, swarm_file, data_file, swarm_path):
    """
    Set up data structure (graph)
    """
    print("Building data structure")
    tim = time.clock()

    if data_file and swarm_file:
        """
        Initilize graph with data from swarm file and data file
        """
        with open(swarm_file, "rU") as swarm_file:
            with open(data_file, "rU") as data_file:
                G = initGraph(swarm_file, data_file, True)
    elif fasta_file:
        """
        Create temporary data from fasta file using swarm.
        """
        # If swarm_path given
        if swarm_path:
            if swarm_path[-1] == "/":
                swarm_path = swarm_path + "swarm"
            swarm = [swarm_path, "-b", "-d", "1"]
        # else check if swarm binary in /usr/bin
        elif os.path.isfile("/usr/bin/swarm"):
            swarm = ["/usr/bin/swarm", "-b", "-d", "1"]
        # else assume swarm file in same folder as this script
        else:  # elif os.path.isfile("./swarm"):
            swarm = ["./swarm", "-b", "-d", "1"]

        with open(fasta_file, "rU") as fasta_file:
            with tempfile.SpooledTemporaryFile() as tmp_swarm_file:
                with tempfile.SpooledTemporaryFile() as tmp_swarm_results:
                    # Run Swarm
                    try:
                        tim2 = time.clock()
                        popen = subprocess.Popen(swarm,
                                                 stderr=tmp_swarm_results,
                                                 stdout=tmp_swarm_file,
                                                 stdin=fasta_file,
                                                 close_fds=True)
                        popen.wait()
                        print("Time for running swarm: ", time.clock() - tim2)
                    except:
                        print("ERROR: SWARM BINARY NOT FOUND OR NOT WORKING")
                        print("Either supply fasta file and Swarm",
                              "binary OR supply data file and swarm file.")
                        sys.exit()

                    tim2 = time.clock()
                    # rewind to the begining of the file
                    tmp_swarm_results.seek(0)
                    # Create data for graph
                    graph_data = [line.strip().split("\t")[1:4]
                                  for line in tmp_swarm_results
                                  if line.startswith("@")]
                    print("Time for extracting data: ", time.clock() - tim2)
                    # rewind to the begining of the file
                    tmp_swarm_file.seek(0)
                    # Initialize graph with obtained data
                    G = initGraph(tmp_swarm_file, graph_data, False)
    else:
        print("ERROR: NO DATA AND SWARM FILES OR FASTA FILE GIVEN")
        sys.exit(0)

    print("Time:", time.clock() - tim)
    print("\nGraph size:", len(G))

    return G


# *************************************************************************** #
#                                                                             #
#                               Pretty printers                               #
#                                                                             #
# *************************************************************************** #

def manualCutter(G, possibleCuts):
    """
    Manual cutting mode
    """
    print("\nENTER MANUAL CUTTING MODE:")
    print("Instructions: To perform cut press y (yes) or enter, "
          "to keep edge press n (no)")
    print("Possible cuts:")
    finalCuts = []
    for edge in possibleCuts:
        # Ignore candidates if both belongs to same root.
        if G[edge[0]].belongingRoot != G[edge[1]].belongingRoot:
            print(G[G[edge[0]].belongingRoot].abundance, "\t\t ",
                  G[G[edge[1]].belongingRoot].abundance)
            print("    \ \t\t /")
            print("    ", G[edge[0]].abundance,
                  "-------", G[edge[1]].abundance)
            minPeak = min(G[G[edge[0]].belongingRoot].abundance,
                          G[G[edge[1]].belongingRoot].abundance)
            maxPeak = max(G[G[edge[0]].belongingRoot].abundance,
                          G[G[edge[1]].belongingRoot].abundance)
            weakPoint = min(G[edge[0]].abundance, G[edge[1]].abundance)
            print("Valley ratio:", float(minPeak) / weakPoint)
            print("Peak ratio  :", float(maxPeak) / minPeak)
            answer = raw_input("Make cut between " + str(G[edge[0]].abundance)
                               + " and " + str(G[edge[1]].abundance) + "? ")
            if not answer == "n" or answer == "no" or answer == "nn":
                finalCuts.append(edge)
    return finalCuts


def prettyPrintCuts(G, finalCuts):
    """
    A bad trial to give a nice output of the cuts.
    """
    print("[Cand for cut]    [Cand abundance]    ", end=" ", file=sys.stderr)
    print("[Closest root]    [Roots abundance]", file=sys.stderr)
    space1 = 18
    space2 = 20
    space3 = 18
    for edge in finalCuts:
        len2 = space1 - len(str([G[node].num for node in edge]))
        len3 = space2 - len(str([G[node].abundance for node in edge])) - 4
        len4 = space3 - len(str([G[node].belongingRoot for node in edge])) - 1
        print(" ", [G[node].num for node in edge], " " * len2,
              [G[node].abundance for node in edge], " " * len3,
              [G[node].belongingRoot for node in edge], " " * len4,
              [G[G[node].belongingRoot].abundance for node in edge],
              file=sys.stderr)
        # print [G[node].name for node in edge]


# *************************************************************************** #
#                                                                             #
#                            Tools for break swarm                            #
#                                                                             #
# *************************************************************************** #

def assignParent(G, threshold):
    """
    Assign biggest neighbour as parent for each node in graph
    """
    print("\nAssigning parents", file=sys.stderr)
    tim = time.clock()

    for node in G:
        # If node is a leaf, set only neighbour as parent if node is small
        if len(node.neighbours) == 1 and node.abundance < threshold:
            node.parent = node.neighbours[0]
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
            # else if no neighbour with degree higher than itself, set as root
            else:
                node.parent = -2  # Has no parent
                node.belongingRoot = node.num

    print("Time:", time.clock() - tim, file=sys.stderr)

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
    print("\nFinding possible cuts...", end="", file=sys.stderr)
    tim = time.clock()

    possibleCuts = []
    for node in G:
        for neighbour in node.neighbours:
            if node.num < neighbour:
                if node.parent != neighbour and \
                   G[neighbour].parent != node.num:
                    possibleCuts.append((node.num, neighbour))
    # print "Possible cuts:\n",possibleCuts
    print("found", len(possibleCuts), "possible cuts", file=sys.stderr)
    # prettyPrintCuts(G,possibleCuts,tim)

    print("Time:", time.clock() - tim, file=sys.stderr)
    return possibleCuts


def findBelongingRoot(G, possibleCuts):
    """
    Find corresponding root for each node in possible cuts
    """
    print("\nFinding belonging roots", file=sys.stderr)
    tim = time.clock()

    for edge in possibleCuts:
        for node in edge:
            if G[node].belongingRoot == -1:  # If -2 then already root
                currentNode = G[node]
                # path = [node]
                while G[currentNode.parent].belongingRoot == -1:
                    # path.append(currentNode.parent)

                    # Skip pair of nodes pointing at each other
                    if G[currentNode.parent].parent == currentNode.num:
                        G[currentNode.parent].belongingRoot = currentNode.num
                        break
                    currentNode = G[currentNode.parent]
                # for n in path:
                #    G[n].belongingRoot = G[currentNode.parent].belongingRoot
                G[node].belongingRoot = G[currentNode.parent].belongingRoot

    print("Time:", time.clock() - tim, file=sys.stderr)


def rewireNode(G, possibleCuts, threshold, root=True):
    """
    If belonging root is below threshold, rewire node from breaking
    point's belonging root or from breaking point. Performs
    breadth-first search (BFS) from starting point to find new
    belonging root. Testing which version is better - rewire from root
    seems most promising.
    """
    print("\nRewiring nodes", file=sys.stderr)
    tim = time.clock()

    for possibleCut in possibleCuts:
        for node in possibleCut:
            if G[G[node].belongingRoot].abundance < threshold:
                newParent = -1
                distance = sys.maxsize
                belongingRootAbundance = 0
                starting_point = G[node].belongingRoot if root else node
                queue = deque([(0, starting_point)])
                while newParent < 0:
                    dist, point = queue.popleft()
                    for neighbour in G[point].neighbours:
                        curAbundance = G[G[neighbour].belongingRoot].abundance
                        if curAbundance > threshold \
                                and curAbundance > belongingRootAbundance \
                                and dist <= distance:
                            newParent = neighbour
                            belongingRootAbundance = curAbundance
                            distance = dist
                        else:
                            queue.append((dist + 1, neighbour))
                # G[node].parent = newParent
                G[node].belongingRoot = G[newParent].belongingRoot
    print("Time:", time.clock() - tim, file=sys.stderr)


def findFinalCuts(G, possibleCuts, threshold, manualCut, parameters):
    """
    Find final cuts, either by manually deciding the cuts or by
    using parameters, or only using the threshold as tiebreaker.
    """
    print("\nFinding final cuts:", file=sys.stderr)
    # Measure time#
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
                # TEST PURPOSE ONLY
                if parameters:
                    weakSpot = min(G[edge[0]].abundance,
                                   G[edge[1]].abundance)
                    biggestRoot = max(G[G[edge[0]].belongingRoot].abundance,
                                      G[G[edge[1]].belongingRoot].abundance)
                    smalletsRoot = min(G[G[edge[0]].belongingRoot].abundance,
                                       G[G[edge[1]].belongingRoot].abundance)
                    if (smalletsRoot / weakSpot > 25
                        and biggestRoot / smalletsRoot < 10) \
                            or smalletsRoot / weakSpot > 50:
                        finalCuts.append(edge)
                else:
                    finalCuts.append(edge)

    # Print Final Cuts - if not too many
    if finalCuts and len(finalCuts) < 100:
        prettyPrintCuts(G, finalCuts)
    print("Number of final cuts:", len(finalCuts), file=sys.stderr)
    print("Time:", time.clock() - tim, "\n", file=sys.stderr)

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
    print("Updating final cuts in data structure:", time.clock() - tim,
          file=sys.stderr)
    # print "Seeds:\n", new_swarm_seeds
    print("Number of possible seeds:", len(new_swarm_seeds), "\n",
          file=sys.stderr)

    return new_swarm_seeds


# *************************************************************************** #
#                                                                             #
#                                 Break swarm                                 #
#                                                                             #
# *************************************************************************** #

def breakSwarm(G, threshold, manualCut, parameters):
    """
    Compute final cuts in the graph.
    """
    # Set THRESHOLD value
    # Default: Ignore roots below 100.
    # if threshold >= 1, then we ignore roots below threshold,
    # else we ignore a percentage of the whole network size.
    threshold = threshold if threshold >= 1 else threshold * len(G)

    # Assign a parent (biggest neighbour) to each node
    assignParent(G, threshold)

    # Find possible cuts
    possibleCuts = findPossibleCuts(G)

    # Find belonging root to each node in possibleCuts
    findBelongingRoot(G, possibleCuts)

    # rewire nodes with small belonging roots
    rewireNode(G, possibleCuts, threshold)

    # find final cuts: manually, parameter, or only by threshold
    finalCuts = findFinalCuts(G, possibleCuts, threshold,
                              manualCut, parameters)

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

    # Update data structure with final cuts
    new_swarm_seeds = updateDataStructure(G, finalCuts)

    return new_swarm_seeds


def findNewSwarms(G, seeds):
    """
    Performing breadth-first search (BFS) from seeds to find new swarms
    """
    print("Performing BFS to discover new swarms...", file=sys.stderr)
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
                new_swarm.append((currentNode, G[currentNode].abundance,
                                  G[currentNode].name))
                for node in G[currentNode].neighbours:
                    if not G[node].checked:
                        queue.append(node)
                        G[node].checked = True

            # Sort amplicons by decreasing abundance and alphabetical order
            new_swarm.sort(key=itemgetter(1, 2), reverse=True)
            new_swarms.append(new_swarm)
    print("Visited nodes:", count, file=sys.stderr)
    print("BFS time:", time.clock() - tim, file=sys.stderr)
    print("Num swarms:", len(new_swarms), file=sys.stderr)
    if count != len(G):
        print("ERROR: DATA LOST AFTER CUTTING.", file=sys.stderr)
    return new_swarms


# *************************************************************************** #
#                                                                             #
#                                     Main                                    #
#                                                                             #
# *************************************************************************** #
def main():
    """
    Main method of the program
    """
    totim = time.clock()

    # Parse command line options
    fasta_file, swarm_file, data_file, threshold, swarm_path, manualCut, \
        parameters = option_parse()

    # Build data structure
    # G = buildGraph(fasta_file, swarm_file, data_file, swarm_path)

    if fasta_file:
        all_amplicons = fasta_parse(fasta_file)

    swarms = swarm_parse(swarm_file)

    # print(swarms)

    all_new_swarms = []
    for swarm in swarms:
        top_amplicon, swarm_mass, swarm_size, top_abundance, amplicons = swarm
        if swarm_size > 2:
            if fasta_file:
                # Run swarm to get the pairwise relationships
                graph_raw_data = run_swarm(swarm_path, all_amplicons,
                                           amplicons, threshold)
                # Build the graph of pairwise relationships
                G = build_graph(graph_raw_data, amplicons, False)
            else:
                graph_raw_data = data_file
                # Build the graph of pairwise relationships
                G = build_graph(graph_raw_data, amplicons, True)

            # Compute cuts and break swarm
            new_swarm_seeds = breakSwarm(G, threshold, manualCut, parameters)

            # Find new swarms
            new_swarms = findNewSwarms(G, new_swarm_seeds)

            # Add to list of all new swarms found
            for swarm in new_swarms:
                all_new_swarms.append(swarm)

    # Output new swarm file
    outputSwarmFile(G, all_new_swarms, swarm_file, fasta_file)

    print("Total time used:", time.clock() - totim, file=sys.stderr)


if __name__ == '__main__':

    main()

    sys.exit(0)


"""
Run example:
python breakOTUs.py -d file.data -s file.swarm
"""

# python swarm_breaker2.py
# -s ./examples/OTU_006_b970fcbdd71ad2a333f702c7ecfe7114.swarm
# -d ./examples/OTU_006_b970fcbdd71ad2a333f702c7ecfe7114.data
# Building data structure
# Network size: 2480
# Time: 0.013227

# Assigning parents
# Time: 0.000986

# Finding possible cuts... found 13 possible cuts
# Time: 0.001195

# Finding belonging roots
# Time: 3.3e-05

# Rewiring nodes
# Time: 5.6e-05

# Finding final cuts:
# [Cand for cut]    [Cand abundance]     [Closest root]    [Roots abundance]
#   [1994, 2171]        [9, 2]            [2169, 2280]       [354, 311]
# Number of final cuts: 1
# Time: 1.4e-05

# Updating final cuts in data structure: 6e-06
# Number of possible seeds: 3

# Performing BFS to discover new swarms...
# Visited nodes: 2480
# BFS time: 0.003767
# Num swarms: 2
# Time used to make output files: 0.003382
# Total time used: 0.023511
