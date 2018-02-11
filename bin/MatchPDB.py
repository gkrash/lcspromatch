#!/usr/bin/env python3

# <Ignore this for now, it's just a placeholder, hasn't been posted anywhere>
# TODO: fill the header section out appropriately.
# GCS-compare for protein structures
# Written by: Gary Krasovic <gkrasovic@gmail.com> in conjunction wth PSH
#

# Imports!
# For squaring and the like.
import math

# For pretty-printing output only.
import pprint
import numpy as np


# Needed for command line stuff
from sys import argv


# Example of what the input file / format looks like:
# We want to skip everything but the backbone...
# format is
# COLUMNS        DATA  TYPE    FIELD        DEFINITION
# -------------------------------------------------------------------------------------
#  1 -  6        Record name   "ATOM  "
#  7 - 11        Integer       serial       Atom  serial number.
# 13 - 16        Atom          name         Atom name.
# 17             Character     altLoc       Alternate location indicator.
# 18 - 20        Residue name  resName      Residue name.
# 22             Character     chainID      Chain identifier.
# 23 - 26        Integer       resSeq       Residue sequence number.
# 27             AChar         iCode        Code for insertion of residues.
# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
# 55 - 60        Real(6.2)     occupancy    Occupancy.
# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
# 77 - 78        LString(2)    element      Element symbol, right-justified.
# 79 - 80        LString(2)    charge       Charge  on the atom.
# Example line:
# ATOM     10  CA  GLU A   2      97.426  76.040  91.645  1.00  3.69      2KAU 374
# Need to really only be concerned with the CA atoms, and the x,y, and Z coords.

# Input: identify the protein we're looking at.
# proteinId1 = '2kau'
# proteinId2 = '1fkx'
# proteinId2 = '2kau'
# inpath = '../pdb/'
# inpath = './'


# Pull our input files from the command line directly..
script, proteinId1, proteinId2 = argv

# This is our error limit, in angstroms beyond which we will not extend further.
# Originally I had 13 angstroms in here, but that's WAY too big.. May be worth stating that if our gap is bigger than 13, we have an error state..

# So I pulled this number from some statistical analysis on the 1AON_A and 1AON_H proteins
# It's the mean difference in intra-CA distance (0.023567944) plus a standard deviation (0.018285768)
# Theory is that if we have a situation where it's a standard deviation above the mean, we're in a weird situation and should gap.
errLim = .04185371 # experimetally derived from the test proteins -


# This is our error function, since we will likely be modifying that some:
def distFunc(x1,y1,z1,x2,y2,z2):
    return(math.pow(math.pow((x1-x2),2) + math.pow((y1-y2),2) + math.pow((z1-z2),2), .5))

# A quick sub to get "next" and "previous" info for a given CA.. called at the end.
def getAtomInfo(protein, number):
    myID = number - 1
    name = protein['name']
    index = protein['serial'][myID]
    residue = protein['residueSeq'][myID]

    p1distanceFromPrevious = 0
    p1distanceFromNext = 0
    p2distanceFromPrevious = 0
    p2distanceFromNext = 0

    if(myID > 0):
        p1distanceFromPrevious = p1pwm[myID][myID-1]
        p2distanceFromPrevious = p2pwm[myID][myID-1]

    if(myID < protein['size'] - 1):
        p1distanceFromNext = p1pwm[myID][myID+1]
        p2distanceFromNext = p2pwm[myID][myID+1]

    print("Protein Name: ",name)
    print("Index: ",index)
    print("Residue Number: ", residue)
    print("Distances from Previous: ", p1distanceFromPrevious, ", ", p2distanceFromPrevious)
    print("Distances from Next: ", p1distanceFromNext, ", ", p2distanceFromNext)


# Sub to print out a list of intra-CA distances, along with a sum of the error for each..
# Maybe this will provide some numberic insight as to what we want to see.
def getDistanceInfo():
    p1name = proteinId1
    p2name = proteinId1

    totDist = 0
    for idx in range(1,protein1['size']):
        resi = protein1['residueSeq'][idx]

        # TODO: Eventually may want to fix this so it works on a "generic" protein so they don't have to be the same length
        distFromPrevResi1 = p1pwm[idx-1][idx]
        distFromPrevResi2 = p2pwm[idx-1][idx]
        # totDist = totDist + distToNextResi
        print(str(resi) + "," + str(distFromPrevResi1) + "," + str(distFromPrevResi2))


# I'm doing this  with separate associatave arrays, because it's simpler for now.
protein1 = {}
protein2 = {}

protein1['name'] = proteinId1
protein1['serial'] = []
protein1['x'] = []
protein1['y'] = []
protein1['z'] = []
protein1['size'] = 0

protein2['name'] = proteinId2
protein2['serial'] = []
protein2['x'] = []
protein2['y'] = []
protein2['z'] = []
protein2['size'] = 0


# Here lies our load code - again, obtuse, but for simplicity.

# infile = inpath + proteinId1 + '.brk'
infile = proteinId1

with open(infile, 'r') as f:

    c = 0 # Counter variable

    for line in f:
        # If the line is an atomic backbone structure line...

        if (line[0:5] == 'ATOM ' and line[12:16] == ' CA '):
            serial = "".join(line[6:11].split())
            x = float("".join(line[30:37].split()))
            y = float("".join(line[38:45].split()))
            z = float("".join(line[46:53].split()))
            protein1['serial'].append(serial)
            protein1['x'].append(x)
            protein1['y'].append(y)
            protein1['z'].append(z)
            c = c + 1 # Increment the counter by 1.
            protein1['size'] = c
f.close()

# infile = inpath + proteinId2 + '.brk'
infile = proteinId2

with open(infile, 'r') as f:

    c = 0 # Counter variable

    for line in f:
        # If the line is an atomic backbone structure line...

        if (line[0:5] == 'ATOM ' and line[12:16] == ' CA '):
            serial = "".join(line[6:11].split())
            x = float("".join(line[30:37].split()))
            y = float("".join(line[38:45].split()))
            z = float("".join(line[46:53].split()))
            protein2['serial'].append(serial)
            protein2['x'].append(x)
            protein2['y'].append(y)
            protein2['z'].append(z)
            c = c + 1 # Increment the counter by 1.
            protein2['size'] = c
f.close()

# At this point our proteins are loaded, and we can continue on to produce our
# matricies

#Protein1PairwiseMatrix is called p1pwm - convention follows for Protein 2
# I'm still being somewhat obtuse with how we're doing this, I depended arrays a lot
# simply because they're fast to address and add to.


p1pwm = []
for row in range(0,protein1['size']):
    p1pwm.append([])
    for col in range(0,protein1['size']):
        x1 = protein1['x'][row]
        y1 = protein1['y'][row]
        z1 = protein1['z'][row]
        x2 = protein1['x'][col]
        y2 = protein1['y'][col]
        z2 = protein1['z'][col]
        p1pwm[row].append(errorFunc(x1,y1,z1,x2,y2,z2))

p2pwm = []
for row in range(0,protein2['size']):
    p2pwm.append([])
    for col in range(0,protein2['size']):
        x1 = protein2['x'][row]
        y1 = protein2['y'][row]
        z1 = protein2['z'][row]
        x2 = protein2['x'][col]
        y2 = protein2['y'][col]
        z2 = protein2['z'][col]
        p2pwm[row].append(errorFunc(x1,y1,z1,x2,y2,z2))

# Now that we have our pairwise matricies, we can produce the greatest common
# subsequence matrices.  The naming convention we will use will be:
# gcsSequence[p1][p2] = [se,qu,en,ce]
# gcsSeqLen[p1][p2] = int(sequenceLength)
# gcsSeqErr[p1][p2] = float(totalSequenceError)
# gcsTlx = float (translation of p2 from x at beginning of chain)
# gcsTly = float (translation of p2 from y at beginning of chain)
# gcsTlz = float (translation of p2 from z at beginning of chain)

gcsSequence = []
gcsSeqLen = []
gcsSeqErr = []
gcsTlx = []
gcsTly = []
gcsTlz = []

for row in range(0,protein1['size']):
    gcsSequence.append([])
    gcsSeqLen.append([])
    gcsSeqErr.append([])
    gcsTlx.append([]) # Translation variables, we're not really using these anymore.
    gcsTly.append([]) # Translation variables, we're not really using these anymore.
    gcsTlz.append([]) # Translation variables, we're not really using these anymore.
    for col in range(0,protein2['size']):

        # Add a sequence array to this cell so we can append to it.
        gcsSequence[row].append([])

        # Now for the cell p1,p2, we need to determine, out of the 4 options,
        # the most locally viable extension, given an error floor of errLim
        p1x = protein1['x'][row]
        p1y = protein1['y'][row]
        p1z = protein1['z'][row]
        p2x = protein2['x'][col]
        p2y = protein2['y'][col]
        p2z = protein2['z'][col]

        # Get Len / Err of the left cell
        if(col > 0): # So we don't go left if we're *on* the left.
            leftErr = gcsSeqErr[row][col-1]
            leftLen = gcsSeqLen[row][col-1]
        else:
            leftErr = float(errLim*2)
            leftLen = 0

        # Get Len / Err of the upper cell
        if(row > 0): # So we don't go left if we're *on* the left.
            upErr = gcsSeqErr[row-1][col]
            upLen = gcsSeqLen[row-1][col]
        else:
            upErr = float(errLim*2)
            upLen = 0

        # Get Len / Err of extending from the diagonal
        if(row > 0 and col > 0):
            diagLen = gcsSeqLen[row-1][col-1] + 1

            # So here, we need to compare the shortest intra-atom distance in both molecules
            # "row" is the index of the new atom in p1, "col" is the index of the new atom in p2.
            # Adding that error to that which was previously calculated - we only need to add in a single new atom from each protein, since we have
            # made reasonably optimal decisions previously.

            # Don't need this either.
            # leastList = []
            # leastLocalErr = errLim * 2 # just so I know this is bigger than anything we'll see.

            # This is the bit of code that calculated ALL of the pairwise distances - I'm commenting it out because now we only want to
            # look at adding to the end of the sequence (which in retrospect makes a lot more sense.)

            # for p1index in range(0,row):
            #     # Get the values..
            #     for p2index in range(0,col):
            #
            #         if(p1pwm[row][p1index] >= errLim or p2pwm[col][p2index] >= errLim):
            #             myErr = errLim * 2 # Our universal sign of "don't go here" if either of  our limits is too big.
            #         else:
            #             myErr = abs(p1pwm[row][p1index] - p2pwm[col][p2index])
            #
            #         # This will get us a single "least" error for each atom pair in our subset
            #         if(myErr < leastLocalErr):
            #             myResult = [myErr,[[row,p1index],[col,p2index]]]
            #             leastList.append(myResult)



            # This line figures out the difference between the current and previous backbones in p1, compared with the current and previous backbones in p2.
            myErr = abs(p1pwm[row][row-1] - p2pwm[col][col-1])
            if( myErr > errLim):
                diagLen = 0 # Simple way of saying "don't do this if the error is too great"
            else:
                diagErr = myErr + gcsSeqErr[row-1][col-1]

            if(localError > errLim):
                diagLen = 0 # simple way of ensuring we don't take the diagonal if it's greater than our error limit.

            diagErr = localError + gcsSeqErr[row-1][col-1]
        else:
            diagErr = float(errLim*2)
            diagLen = 0

        # What's it look like if we start again with a new chain
        selfLen = 1
        selfErr = 0

        # Now we need to build out a list of the options so we can max() it

        # print(leftErr, upErr, diagErr, selfErr)
        extension = [[leftLen,leftErr,'left'],[upLen,upErr,'up'],[diagLen,diagErr,'diag'],[selfLen,selfErr,'self']]

        # Now we can do this fanciness to get max length FIRST with min error for that max length.
        # TODO: This is really a place that we can look to make some improvements
        foundLen,foundErr,direction = max(extension, key=lambda x: (x[0],-x[1]))

        gcsSeqLen[row].append(foundLen)
        gcsSeqErr[row].append(foundErr)


        # print(direction,row,col, foundErr)
        # Now we just need to populate the other arrays.
        if(direction == 'left'):
            gcsSequence[row][col] = gcsSequence[row][col - 1]
            # gcsTlx[row].append(gcsTlx[row][col - 1])
            # gcsTly[row].append(gcsTly[row][col - 1])
            # gcsTlz[row].append(gcsTlz[row][col - 1])
        elif(direction == 'up'):
            gcsSequence[row][col] = gcsSequence[row - 1][col]
            # gcsTlx[row].append(gcsTlx[row - 1][col])
            # gcsTly[row].append(gcsTly[row - 1][col])
            # gcsTlz[row].append(gcsTlz[row - 1][col])
        elif(direction == 'diag'):
            gcsSequence[row][col] = []
            gcsSequence[row][col] += gcsSequence[row-1][col-1]

            # print("Copied from Diagnonal " + str(gcsSequence[row][col]))
            gcsSequence[row][col].append([row,col])

            # print("Final value "+ str(gcsSequence[row][col]))
            # gcsTlx[row].append(gcsTlx[row - 1][col])
            # gcsTly[row].append(gcsTly[row - 1][col])
            # gcsTlz[row].append(gcsTlz[row - 1][col])
        elif(direction == 'self'):
            gcsSequence[row][col] = [[row,col]]

            # This is a funny one, since it's a new chain, we're re-centering, and storing what we're translating by so we can use it later.
            # # gcsSeqErr[row][column] = 0 - may consider bringing this back later - I want to leave something here now though to demerit bad matches.
            # gcsTlx[row].append(p1x - p2x)
            # gcsTly[row].append(p1y - p2y)
            # gcsTlz[row].append(p1z - p2z)





# print(np.matrix(gcsSeqLen))
# print(np.matrix(gcsSeqErr))

longestChain = str(len(gcsSequence[protein1['size']-1][protein2['size']-1]))

# Error calculations
e = np.array(gcsSeqErr)
e.flatten()
errorMean = str(np.mean(e))
errorStdDev = str(np.std(e))

# print("Protein 1 backbone length: " + str(protein1['size']))
# print("Protein 2 backbone length " + str(protein2['size']))
# print("Longest common set of backbones: " + str(len(gcsSequence[protein1['size']-1][protein2['size']-1])))
# pprint.pprint(gcsSeqErr)


print(longestChain + "," + errorMean)
