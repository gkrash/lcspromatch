#!/usr/bin/env python3

# <Ignore this for now, it's just a placeholder, hasn't been posted anywhere>
# TODO: fill the header section out appropriately.
# GCS-compare for protein structures
# Written by: Gary Krasovic <gkrasovic@gmail.com> in conjunction wth PSH
#


# Imports!
# For squaring and the like.
import math

# Needed for command line stuff
from sys import argv
import inspect, os


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

# Defining our local tmpdir for output if we need it.
tmpdir = str(os.path.dirname(os.path.dirname(os.path.abspath(inspect.stack()[0][1]))))+"/tmp"

# Pull our input files from the command line directly..
script, proteinId1, proteinId2 = argv

# This is our error limit, in angstroms beyond which we will not extend further.
# Originally I had 13 angstroms in here, but that's WAY too big.. May be worth stating that if our gap is bigger than 13, we have an error state..

# So I pulled this number from some statistical analysis on the 1AON_A and 1AON_H proteins
# It's the mean difference in intra-CA distance (0.023567944) plus a standard deviation (0.018285768)
# Theory is that if we have a situation where it's a standard deviation above the mean, we're in a weird situation and should gap.
# errLim = .04185371 # This was a broken way of doing this..
errLim = 9999 # garbage value for now.

# This defines a local distance for our distance calculations, providing a local limit
localDistCap = 99999 # Eventually we want to flip this back to 13 (angstroms), but not for now.


# This is the function that measures distance between 2 atoms
def distFunc(x1,y1,z1,x2,y2,z2):
    return(math.pow(math.pow((x1-x2),2) + math.pow((y1-y2),2) + math.pow((z1-z2),2), .5))

# I'm doing this  with separate associatave arrays, because it's simpler for now.
protein1 = {}
protein2 = {}

protein1['name'] = os.path.splitext(os.path.basename(proteinId1))[0]
protein1['serial'] = []
protein1['residueSeq'] = []
protein1['x'] = []
protein1['y'] = []
protein1['z'] = []
protein1['size'] = 0
protein1['source'] = [] # we're going to pull the whole thing into an object too

protein2['name'] = os.path.splitext(os.path.basename(proteinId2))[0]
protein2['serial'] = []
protein2['residueSeq'] = []
protein2['x'] = []
protein2['y'] = []
protein2['z'] = []
protein2['size'] = 0
protein2['source'] = [] # we're going to pull the whole thing into an object too


# Here lies our load code - again, obtuse, but for simplicity.

# infile = inpath + proteinId1 + '.brk'
infile = proteinId1

with open(infile, 'r') as f:

    c = 0 # Counter variable

    for line in f:
        # If the line is an atomic backbone structure line...

        if (line[0:5] == 'ATOM ' and line[12:16] == ' CA '):
            protein1['source'].append(line.rstrip())
            serial = "".join(line[6:11].split())
            residueSeq = "".join(line[22:26].split())
            x = float("".join(line[30:37].split()))
            y = float("".join(line[38:45].split()))
            z = float("".join(line[46:53].split()))
            protein1['serial'].append(serial)
            protein1['residueSeq'].append(residueSeq)
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
            protein2['source'].append(line.rstrip())
            serial = "".join(line[6:11].split())
            residueSeq = "".join(line[22:26].split())
            x = float("".join(line[30:37].split()))
            y = float("".join(line[38:45].split()))
            z = float("".join(line[46:53].split()))
            protein2['serial'].append(serial)
            protein2['residueSeq'].append(residueSeq)
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
        # Determine the distance between the two C-alphas (backbones)
        p1pwm[row].append(distFunc(x1,y1,z1,x2,y2,z2))

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
        # Determine the distance between the two C-alphas (backbones)
        p2pwm[row].append(distFunc(x1,y1,z1,x2,y2,z2))

# Now that we have our pairwise matricies, we can produce the greatest common
# subsequence matrices.  The naming convention we will use will be:
# gcsSequence[p1][p2] = [se,qu,en,ce]
# gcsSeqLen[p1][p2] = int(sequenceLength)
# gcsSeqErr[p1][p2] = float(totalSequenceError)

gcsSequence = []
gcsSeqLen = []
gcsSeqErr = []

for row in range(0,protein1['size']):
    gcsSequence.append([])
    gcsSeqLen.append([])
    gcsSeqErr.append([])

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
            leftErr = 0
            leftLen = -1

        # Get Len / Err of the upper cell
        if(row > 0): # So we don't go left if we're *on* the left.
            upErr = gcsSeqErr[row-1][col]
            upLen = gcsSeqLen[row-1][col]
        else:
            upErr = 0
            upLen = -1

        # Get Len / Err of extending from the diagonal
        if(row > 0 and col > 0):
            diagLen = gcsSeqLen[row-1][col-1] + 1

            # So here, we need to compare the shortest intra-atom distance in both molecules
            # "row" is the index of the new atom in p1, "col" is the index of the new atom in p2.
            # Adding that error to that which was previously calculated - we only need to add in a single new atom from each protein, since we have
            # made reasonably optimal decisions previously.

            # This line figures out the difference between the current and previous backbones in p1, compared with the current and previous backbones in p2.
            myErr = abs(p1pwm[row][row-1] - p2pwm[col][col-1])

            if( myErr > errLim):
                diagLen = 0 # Simple way of saying "don't do this if the error is too great"
            else:
                diagErr = myErr + gcsSeqErr[row-1][col-1]

        # slightly confusing, but.. do this if we're in the first row or col..
        else:
            diagErr = 0
            diagLen = -1

        # We also need to check all of these against the current cell as it stands.
        if(row > 0 and col > 0):
            selfLen = 1 # if we pick this, we've got at least a length of 1
            selfErr = abs(p1pwm[row][row-1] - p2pwm[col][col-1])

        # if we're in the first row or col, don't do anything at all..
        else:
            selfErr = 0
            selfLen = 0

        # Now we need to build out a list of the options so we can max() it

        # print(leftErr, upErr, diagErr, selfErr)
        extension = [[selfLen,selfErr,'self'],[leftLen,leftErr,'left'],[upLen,upErr,'up'],[diagLen,diagErr,'diag']]

        # Now we can do this fanciness to get max length FIRST with min error for that max length.
        # TODO: This is really a place that we can look to make some improvements
        foundLen,foundErr,direction = max(extension, key=lambda x: (x[0],-x[1]))

        gcsSeqLen[row].append(foundLen)
        gcsSeqErr[row].append(foundErr)

        # print(direction,row,col, foundErr)
        # Now we just need to populate the other arrays.
        if(direction == 'left'):
            gcsSequence[row][col] = []
            gcsSequence[row][col] += gcsSequence[row][col - 1]

        elif(direction == 'up'):
            gcsSequence[row][col] = []
            gcsSequence[row][col] += gcsSequence[row - 1][col]

        elif(direction == 'diag'):
            gcsSequence[row][col] = []
            gcsSequence[row][col] += gcsSequence[row-1][col-1]
            gcsSequence[row][col].append([row,col])

        elif(direction == 'self'):
            gcsSequence[row][col] = [[row,col]]

    # if(row==0): print(gcsSequence[row])





# output a selection of atoms in PDB format, from each structure so that we can import to PyMOL for analysis
# We'll also pull the coordinates for the matched proteins in both, so that we can do our RMSD calculation on them.
# This will output in the 'tmp' folder in the root of the project.

outfile1 = tmpdir + os.sep + protein1['name'] + "-aligned.pdb"
outfile2 = tmpdir + os.sep + protein2['name'] + "-aligned.pdb"
print(outfile1)

with open(outfile1, 'w') as f1, open(outfile2, 'w') as f2:
    for idx in range(0,len(gcsSequence[protein1['size']-1][protein2['size']-1])):
        p1idx,p2idx = gcsSequence[protein1['size']-1][protein2['size']-1][idx]
        f1.write(protein1['source'][p1idx]+'\n')
        f2.write(protein2['source'][p2idx]+'\n')
f1.close()
f2.close()


    # print(protein1['residueSeq'][p1idx], protein1['residueSeq'][p2idx])


# Error calculations
e = np.array(gcsSeqErr)
e.flatten()
errorMean = str(np.mean(e))
errorStdDev = str(np.std(e))

print("Protein 1 backbone length: " + str(protein1['size']))
print("Protein 2 backbone length " + str(protein2['size']))
# print("Longest common set of backbones: " + str(len(gcsSequence[protein1['size']-1][protein2['size']-1])))
# pprint.pprint(gcsSeqErr)


# Get info about a particular CA Atom...
#getAtomInfo(protein1,267)

# Print out distance info
# getDistanceInfo()
print(chainSize + "," + errorMean)
