# lcspromatch
Matching Proteins with a Longest Common Subsequence Approach

## Things you should know before running this:

## The DALILite build WILL LIKELY NOT WORK
It was installed from the original comparison work we were doing, and since it's
built on a bunch of tightly integrated code, requires a bunch of Make / dependency
work to be installed properly.  If you want to use DaliLite, I suggest you go grab it from their
homepage, and work with it.  I made no modifications to it, just wrote some surround scripts
to pull in data, and parse output.

## Install instructions for the rest:
### Dependencies:
* These scripts depend on Python3 being in your path
* These scripts also depend on `/usr/bin/env` - the first choice here is where we grab our python3 interpreter from.
* These scripts have some required includes.  I handle these with pipenv.

### Installation
1. Clone the git repo
1. Install Python 3.6 if you haven't got it already, make sure you can execute it with `python3`
1. Install pip if you haven't already
1. Install pipenv (to your user dir) if you haven't already `pip3 install --user pipenv`
1. Next, set up the environment with pipenv, cd to the project dir and execute `pipenv --python 3.6`
1. Fire up a shell in this environment with `pipenv shell` and start running things
**note:** if you DO want to install some additional packages, hop out of that shell, and install them with `pipenv install <package>` from the repository root.  This will ensure they get added to the Pipfile.

That should do it, all that needs to happen at this point is to try out the `bin/MatchPDB.py` script with your 2 favorite PDB files.
