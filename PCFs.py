#! /usr/bin/env python

import sys
import re
import subprocess
import os
import shutil
import itertools

"""
Pipeline for generating Phylogeographic Concordance Factors. Input is 
posterior distribution of trees for each species (from *BEAST or SNAPP),
and the script will clean the files, calculate counts for each unique 
topology for each species (mbsum), and call BUCKy to calculate a concordance 
tree and concordance factors. This returns a table with all permutations of
taxa from K = 2 to K = N taxa and their corresponding average PCF value.
A folder with all BUCKy runs is also saved and given a unique model number.
See README.md for specifics of requirements for running the script.

Author: Jordan Satler
Date: 10 April 2015
Version: 1
"""

def clean(file): 
    """removes population parameter information from tree files"""
    Openf = open(file, 'r')
    OutfileName = Openf.name[:-6]
    Out = []
    Trees = 0
    for line in Openf:
        line = line.strip()
        if "STATE" not in line:
            Out.append(line)
        elif "STATE" in line:
            r = re.sub('\[(.*?)\]', '', line)
            Out.append(r)
            Trees += 1
    Openf.close()
    return Out, OutfileName, Trees

def OGnum(postDist, filename):
    """add two outgroups to tree data. necessary for calculating PCFs."""
    MB_file = filename + '_mbsumReady.txt'
    OpenW = open(MB_file, 'w')
    Pops = []
    Traits = {}
    Added = ''

    for line in postDist:
        
        if 'STATE' not in line:
                   
            if re.search('^\d+', line):
                line_a = line.strip(',').split()
                line_b = line.strip(',')
                
                #Turn this into an integer so OG can be double digits.
                Pops.append(int(line_a[0]))
                Added += line_b + ',' + '\n'
                Traits[line_a[0]] = line_a[1]
                
            elif line == ';' and len(Pops) > 0:
                to_add = "%d OG1,\n%d OG2\n;\n" % (int(max(Pops)) + 1, int(max(Pops)) + 2)
                Added += to_add
                
            else:
                Added += line + '\n'

        elif 'STATE' in line:
            s = re.sub(r'(\(.*\))', r'((\1,' + str(OG1) + '),' + str(OG2) + ')', line)
            Added += s + '\n'
        
        if Pops:
            OG1 = int(max(Pops)) + 1
            OG2 = OG1 + 1
    
    Traits[str(int(max(Pops)) + 1)] = 'OG1'
    Traits[str(OG1 + 1)] = 'OG2'
                       
    OpenW.write(Added)
    OpenW.close()
    return MB_file, Traits

def mbsum(trees, PostSum):
    """calls mbsum for unique topology counts. discards first 10% of trees 
       as burnin, but this can be set by the user below. moves files into 
       their respective folders"""
    OpenMB = open(trees, 'r')
    Name = OpenMB.name
    Out = Name[:-4] + "_mbsum_Results.txt"
    
    #Value can be changed to reflect desired burn-in
    Post = PostSum * 0.1
    
    call = subprocess.call(["mbsum", "-n",  str(int(Post)), "-o",  str(Out), str(Name)])
    OpenMB.close()
    
    #Move mbsum input files
    if os.path.exists("./mbsum/mbsum_in/"):
        shutil.move("./" + str(Name), "./mbsum/mbsum_in/" + str(Name))
    else:
        os.makedirs("./mbsum/mbsum_in/")
        shutil.move("./" + str(Name), "./mbsum/mbsum_in/" + str(Name))
    
    #Move mbsum output files
    if os.path.exists("./mbsum/mbsum_out/"):
        shutil.move("./" + str(Out), "./mbsum/mbsum_out/" + str(Out))
    else:
        os.makedirs("./mbsum/mbsum_out")
        shutil.move("./" + str(Out), "./mbsum/mbsum_out/" + str(Out))

def bucky():
    """This function will take in summaries of tree topology distributions
       and create a concordance tree with concordance factors."""
    path = "./mbsum/mbsum_out/"
    trees = [filename for filename in os.listdir(path)]
    
    #This section returns community tree for all permutations
    Data_sets = combos(trees)
    Num = 1
    for j in Data_sets[0]:
        cmd = ["bucky", "--use-independence-prior", "-n", "100000", "-o", "PCF"]
        
        for l in range(len(j)):
            full = path + j[l]
            cmd.append(full)
        runBucky = subprocess.call(cmd)    
    
        if not os.path.exists("./bucky_results/"):
            os.mkdir("./bucky_results")
        if not os.path.exists("./bucky_results/bucky_" + str(Num) + "/"):
            os.makedirs("./bucky_results/bucky_" + str(Num))
        
        Out_bucky = [file for file in os.listdir("./")]
        if not os.path.exists("./bucky_results/bucky_" + str(Num) + "/"):
            os.makedirs("./bucky_results/bucky_" + str(Num))
        
        for i in Out_bucky:
            if 'PCF' in i and i != 'PCF.concordance':
                shutil.move(str(i), "./bucky_results/bucky_" + str(Num) + "/" + str(i))
            
            elif 'PCF' in i and i == 'PCF.concordance':    
                buildPCF = concordance_tree(Traits)
                shutil.move(str(i), "./bucky_results/bucky_" + str(Num) + "/" + str(i))
                
                nodalSup = calculate()
                       
                taxa = []
                for tax in range(len(j)):
                    taxonName = j[tax]
                    taxa.append(taxonName[:-29])
                models = ', '.join(taxa)
                    
                with open("All_Combinations.txt", 'a') as All:
                    if Num == 1:
                        Header = "Model\tK\tAverage\tTaxa"
                        Out = "%d\t%d\t%.4f\t%s" % (Num, len(j), nodalSup, models)
                        All.write(Header + '\n' + Out + '\n')
                    if Num > 1:
                        Out = "%d\t%d\t%.4f\t%s" % (Num, len(j), nodalSup, models)
                        All.write(Out + '\n')

                shutil.move("PCF_Tree.tre", "./bucky_results/bucky_" + str(Num) + "/PCF_Tree.tre")
        Num += 1
    return trees
    
def combos(trees):
    """This returns all possible permutations of taxa, from K of 2 
       to K of length(taxa), as a list."""
    comb = []
    for i in range(2, len(trees) + 1):
        for subset in itertools.combinations(trees, i):
            if subset not in comb:
                comb.append(subset)
    return comb, len(trees)

def concordance_tree(Traits):
    """Takes Bucky output and returns a modified concordance tree with 
        concordance factors"""
    Infile = open("./PCF.concordance", 'r')
    
    lines = []
    tmp = ''
    for line in Infile:
        line = line.strip()
        if tmp:
            lines.append(line)
            tmp = ''
        if line == 'Primary Concordance Tree with Sample Concordance Factors:':
            tmp = line
    
    Tree = ''
    for i in lines:
        tmp = ''
        for j in range(len(i)):
            #Add in OTU names
            if i[j] == ':' and tmp.isdigit() == True:
                add = Traits[tmp] + ":"
                Tree += add
            
            elif i[j].isdigit() == True:
                tmp += i[j]
                continue
            
            else:
                tmp += i[j]
                Tree += tmp
            
            tmp = ''
    
    with open('PCF_Tree.tre', 'w') as PCF:
        PCF.write(Tree[1:len(Tree) - 21] + ');')    
        
    Infile.close()

def calculate():
    """This will calculate the average nodal support"""
    with open("PCF_Tree.tre", 'r') as tree:
        for line in tree:
            line = line.strip()
            
            x = re.findall("\):\d.\d+", line)

            Total = 0
            for i in range(len(x) - 1):
                Total += float(x[i][2:])
            return Total / (len(x) - 1)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "python PCFs.py Input*"
        sys.exit()
    Filelist = sys.argv[1:]
    for file in Filelist:
        Out, OutfileName, Trees = clean(file)
        MB_file, Traits = OGnum(Out, OutfileName)
        mb_out = mbsum(str(MB_file), Trees)
        
        if Filelist.index(file) == len(Filelist) - 1:  
            PCFs = bucky()