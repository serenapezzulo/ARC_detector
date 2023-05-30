ddsim --compactFile ./compact/arc_full_v0.xml --runType run --part.userParticleHandler='' --macroFile myscripts/overlap.mac > overlapDump.txt
grep -A 9 'G4Exception-START' overlapDump.txt
