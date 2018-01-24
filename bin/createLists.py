#!/usr/bin/python3

import os
from random import shuffle

# Script to create 3 separate lists of proteins, and pairwise comparisons within them
# for use in other programs that we're building.
#
# Output is in the form of 3 lists, stored in the root/lists directory, in the form of <path>,<path>


# Directory containing our structure files.
# inpath = '../pdb/'
inpath = '../../pdb/pdbstyle-2.06'

# Where we'll stick these lists when we're done.
outpath = '../'

# How many lists do we want?
numLists = 3

# Do we want these lists to have even numbers of proteins? (1 = true, 0=false)
evenLists = 1

# OPTIONAL - impose a hard cap on the length of our lists
listLen = 0

proteinList = []
eachList = 0
finalLists = []

# Read in our list of proteins and randomize the order.
for dirpath, dirs, files in os.walk(inpath):
    for file in files:
        proteinList.append(os.path.join(dirpath, file))

shuffle(proteinList)

if(listLen > 0):
    eachList = listLen
else:
    # Figure out how long our lists need to be..
    totalNum = len(proteinList)
    eachList = totalNum // numLists # (so we make sure we've got an even number in each one.)

if((evenLists == 1) and (eachList % 2 == 1)):
    eachList -= 1


for l in range(0, numLists):
    myListName = outpath + 'pList' + str(l)
    myFH = open(myListName, 'w')
    myList = [proteinList[i:i + eachList] for i in range(0, len(proteinList), eachList)][l]
    for item in myList:
        myFH.write("%s\n" % item)
#
# for l in range(0, numLists):
#     print(len([proteinList[i:i + eachList] for i in range(0, len(proteinList), eachList)][l]))




print("Produced " + str(numLists) + " lists of length " + str(eachList) + " from a source of " + str(len(proteinList)) + " proteins.")
