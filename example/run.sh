#!/bin/sh
runquick.cuda ACE_ALA10_NME.in 
runquick.cuda ACE_GLY6_NME.in
runquick.cuda taxol.in
runquick ACE_ALA10_NME.in 
runquick ACE_GLY6_NME.in
runquick taxol.in

