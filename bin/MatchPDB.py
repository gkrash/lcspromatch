#!/usr/bin/env python3

# <Ignore this for now, it's just a placeholder, hasn't been posted anywhere>
# TODO: fill the header section out appropriately.
# GCS-compare for protein structures
# Written by: Gary Krasovic <gkrasovic at gmail.com> in conjunction wth PSH
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

# Converted the whole protein init / import process to a sub, hopefully this will make it easier to read.
def initProtein(proteinId):

    myProtein = {} # Initialize the map

    # Init the variables stored in the map
    myProtein['name'] = os.path.splitext(os.path.basename(proteinId))[0] # The name of the protein dervied from the filename
    myProtein['serial'] = [] # Array of serial numbers of c-alphas found
    myProtein['residueSeq'] = [] # Array of residue numbers of c-alphas found
    myProtein['x'] = [] # Array of x-coords of c-alphas found
    myProtein['y'] = [] # Array of y-coords of c-alphas found
    myProtein['z'] = [] # Array of z-coordsd of c-alphas found
    myProtein['size'] = 0 # Size of the protein (in number of c-alphas found)
    myProtein['source'] = [] # Source of all c-alphas found
    myProtein['distMat'] = [] # The internal protein pairwise distance matrix.. we only initialize it here.

    myProtein['matchedSubsequence'] = [] # The subsequence matched within the protein, we'll populate this as we go

    with open(proteinId, 'r') as f:

        c = 0 # Counter variable

        for line in f:
            # If the line is an atomic backbone structure line...

            if (line[0:5] == 'ATOM ' and line[12:16] == ' CA '):
                # Read the variables in
                myProtein['source'].append(line.rstrip())
                serial = "".join(line[6:11].split())
                residueSeq = "".join(line[22:26].split())
                x = float("".join(line[30:37].split()))
                y = float("".join(line[38:45].split()))
                z = float("".join(line[46:53].split()))

                # And append them to the protein
                myProtein['serial'].append(serial)
                myProtein['residueSeq'].append(residueSeq)
                myProtein['x'].append(x)
                myProtein['y'].append(y)
                myProtein['z'].append(z)
                c = c + 1 # Increment the counter by 1 - size will be len+1
                myProtein['size'] = c
    f.close()

    return myProtein

# Sub to calcuate error threshold, for easy editing.
def calcErrorThreshold(seqLen):
    return (((seqLen + 1) * (seqLen / 2)) * (1.1 * 1.1))

# Sub to calculate the pairwise distance matrix for a given protein.
def calcPairwiseDistanceMatrix(protein):

    for row in range(0,len(protein['serial'])):
        protein['distMat'].append([])
        for col in range(0,len(protein['serial'])):
            x1 = protein['x'][row]
            y1 = protein['y'][row]
            z1 = protein['z'][row]
            x2 = protein['x'][col]
            y2 = protein['y'][col]
            z2 = protein['z'][col]
            # Determine the distance between the two C-alphas (backbones)
            protein['distMat'][row].append(distFunc(x1,y1,z1,x2,y2,z2))

# Read in the two proteins that we have.
p1 = initProtein(proteinId1)
p2 = initProtein(proteinId2)

# Now we need to calculate the internal distance matrix within protein 1 and protein 2
calcPairwiseDistanceMatrix(p1)
calcPairwiseDistanceMatrix(p2)



gcss = [] # gcss is Greatest Common SubSequence
gcssErr = []

# Every c-alpha in p1 is a 'row' in our gcss table
for row in range(0,len(p1['serial'])):
    print("CA " + str(row+1) + " of " + str(len(p1['serial'])+1))
    gcss.append([])
    gcssErr.append([])

    # Every c-alpha in p2 is a 'col' in our gcss table
    for col in range(0,len(p2['serial'])):

        # Create the column (and thereby the cell)
        gcss[row].append([])
        gcssErr[row].append([])

        # just to get this out of the way, we get the up and side match-lengths

        # Init to zero so we don't choose them if we're in a 0 row/col
        uLen = 0
        sLen = 0

        if(row > 0):
            uLen = len(gcss[row-1][col])
        if(col > 0):
            sLen = len(gcss[row][col-1])

        # Now, for the diagonal

        dErr = 0 # Error for the match stored in the diagnoal
        dLen = 0 # Length of the match stored in the current diagnoal
        dSeq = [] # Sequence from the diagonal
        nErr = 0 # The new error from adding in in next pair of carbon-alphas from p1 and p2

        # Make sure we're in a cell that *has* a diagonal to extend from
        if((col > 0) and (row > 0)):

            # gcss[row-1][col-1] is the matched subsequence from the upper right cell, as a list of pairs
            # dSeq is the actual match sequence stored in the diagonal
            dSeq = gcss[row-1][col-1]

            dLen = len(dSeq)

            # gcssErr[row-1][col-1] is the currently matched subsequence error from the diagonal
            dErr = gcssErr[row-1][col-1]

            # Calculate our error threshold
            eThresh = calcErrorThreshold(dLen)

            # This calculates nErr
            for idx in range(0,dLen):

                p1Atom,p2Atom = dSeq[idx] # Pull the internal indexes for each atom within our aligned structures.

                # calculate the distance between the current atoms we want to add and the one we're targeting from the existing structure.
                p1AtomDist = p1['distMat'][row][p1Atom]
                p2AtomDist = p2['distMat'][col][p2Atom]

                # Only want to add to the error if we're lower than the localDistCap
                if(p1AtomDist <= localDistCap and p2AtomDist <= localDistCap):
                    nErr += math.pow(abs(p1AtomDist - p2AtomDist), 2) # Finally, compare the distances between the two atoms and add the square of that in.
                else:
                    continue



        # Calculate the error threshold..
        errThresh = calcErrorThreshold(dLen)

        # If the error is less than our threshold, extend the existing alignment from the diagnoal, adding in the new row,col
        # Because of the way that dErr and nErr are designed, if this is a 0 row or col,
        if((dErr + nErr) <= errThresh):
            gcssErr[row][col] = dErr + nErr
            gcss[row][col] = dSeq[:] + [[row, col]] # Python normally operates on lists by reference, the colon makes it copy..

        # Or else, if the side alignment length is longer than the one above, copy it  (and it's error) into the current cell.
        elif(uLen < sLen):
            gcssErr[row][col] = gcssErr[row][col-1]
            gcss[row][col] = gcss[row][col-1][:] # Python normally operates on lists by reference, the colon makes it copy..

        # Lastly, if nothing else, copy from above.
        else:
            gcssErr[row][col] = gcssErr[row-1][col]
            gcss[row][col] = gcss[row-1][col][:] # Python normally operates on lists by reference, the colon makes it copy..


# output a selection of atoms in PDB format, from each structure so that we can import to PyMOL for analysis
# We'll also pull the coordinates for the matched proteins in both, so that we can do our RMSD calculation on them.
# This will output in the 'tmp' folder in the root of the project.

outfile1 = tmpdir + os.sep + p1['name'] + "-aligned.pdb"
outfile2 = tmpdir + os.sep + p2['name'] + "-aligned.pdb"

with open(outfile1, 'w') as f1, open(outfile2, 'w') as f2:
    for idx in range(0,len(gcss[p1['size']-1][p2['size']-1])):
        p1idx,p2idx = gcss[p1['size']-1][p2['size']-1][idx]
        f1.write(p1['source'][p1idx]+'\n')
        f2.write(p2['source'][p2idx]+'\n')
f1.close()
f2.close()

# Calculate and print the max length chain size
chainSize = str(len(gcss[-1][-1])) # a -1 index means "last" in python
print("Alignment Size: ", chainSize)
print("Alignments saved as:")
print(outfile1)
print(outfile2)
