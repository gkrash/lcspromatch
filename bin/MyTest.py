#!/usr/bin/python3

import os, sys, shutil, subprocess

from sys import argv


# Script to test comparing 2 proteins with our own technique.

# Execution is in the form of "MyTest.py <protein1> <protein2>"
# Input files are expected to follow the PDB95 standard from <TODO: Insert citation to site I downloaded the source from>

# Output is in the form of: longest-chain, error-mean, true-false-match

# Where:

# the longest-chain and error-mean come from my script above
# The heierarchy match is at the following level:

# The heierarchy match between the two is true or false depending on the variable below.

# Match Level corresponds to:
# https://scop.berkeley.edu/help/ver=1.71
# where:
#
# 0 = Root, all should match here (useless)
# 1 = Class
# 2 = Fold
# 3 = Superfamily
# 4 = Family

matchLevel = 3

# myPDB Matcher Location
myCMD = "/home/krash/git/gcs-comp/bin/MatchPDB.py"

# Temp space
tmp = "/tmp"


## Start program
# Return Variables:
meanerr = float(0)
lchain = 0
isMatch = 0

# Command line input
script, protein1File, protein2File = argv

# what to match on for the heierarchy from SCOP:
matchOn = "SCOPe-sccs:"
nameLine = "SCOPe-sid:"

protein1Name = ""
protein1Categorization = ""
with open(protein1File, 'r') as f:
    for line in f:
        if(len(line.split()) >= 4 and line.split()[3] == matchOn):
            protein1CatList = line.split()[4].split('.')
        elif(len(line.split()) >= 4 and line.split()[3] == nameLine):
            protein1Name = line.split()[4]

protein2Name = ""
protein2Categorization = ""
with open(protein2File, 'r') as f:
    for line in f:
        if(len(line.split()) >= 4 and line.split()[3] == matchOn):
            protein2CatList = line.split()[4].split('.')
        elif(len(line.split()) >= 4 and line.split()[3] == nameLine):
            protein2Name = line.split()[4]

# Compare the two categorizations first and compare to the order we want.
for i in range(0,matchLevel):
    # print(protein1CatList[i],protein2CatList[i])
    if(protein1CatList[i] != protein2CatList[i]):
        isMatch = 0
        break
    else:
        isMatch = 1

# Run the our program..
proc = subprocess.Popen([myCMD,protein1File,protein2File],stdout=subprocess.PIPE)

while True:
  line = proc.stdout.readline().decode()
  if line != '':
    #the real code does filtering here
    print(protein1Name + "," + protein2Name + "," + line.rstrip() + "," + str(isMatch))
  else:
    break
