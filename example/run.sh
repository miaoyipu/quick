#!/bin/sh

mpirun -n 2 quick ACE_ALA10_NME.in
cp ACE_ALA10_NME.out ACE_ALA10_NME_divcon_n2.out

mpirun -n 1 quick ACE_ALA10_NME.in
cp ACE_ALA10_NME.out ACE_ALA10_NME_divcon_n1.out

mpirun -n 4 quick ACE_ALA10_NME.in
cp ACE_ALA10_NME.out ACE_ALA10_NME_divcon_n4.out

mpirun -n 8 quick ACE_ALA10_NME.in
cp ACE_ALA10_NME.out ACE_ALA10_NME_divcon_n8.out

