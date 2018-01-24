#!/usr/bin/python3

import os, sys, shutil, subprocess

from sys import argv


# Wrapper script to run tests across lists of protein chains.


# Execution is in the form of "wrapxxx listFile"

# matchProtein is the protein to test against.
# listFile must contain a list of full paths to protein files.


# Command line input
script, matchProtein, listFile = argv

myCMD = '/home/krash/git/gcs-comp/bin/testcat.py'

myList = []


with open(listFile, 'r') as f:

    c = 0 # Counter variable
    for line in f:
        myList.append(line)
f.close()


for p2 in range(0,len(myList)):
    if(myList[p2].rstrip() == matchProtein):
        next
    else:
        # Run the our program..
        proc = subprocess.Popen([myCMD,matchProtein,myList[p2].rstrip()],stdout=subprocess.PIPE)
        while True:
          line = proc.stdout.readline().decode()
          if line != '':
            #the real code does filtering here
            print(str(line.rstrip()))
          else:
            break
