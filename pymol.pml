fetch 1AON_H 1AON_A
load tmp/1AON_A-aligned.pdb
load tmp/1AON_H-aligned.pdb

align 1AON_A, 1AON_H
align 1AON_H-aligned, 1AON_H
align 1AON_A-aligned, 1AON_A

hide all

show sticks, backbone

show spheres, model 1AON_A-aligned and name CA
show spheres, model 1AON_H-aligned and name CA

alter elem c, vdw=.5

color lightorange, model 1AON_A

color palegreen, model 1AON_H

color magenta, model 1AON_A-aligned and name CA
color cyan, model 1AON_H-aligned and name CA




rebuild
