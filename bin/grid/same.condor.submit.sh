######################################################################
# 				MixedEvent
######################################################################


Rank = ( (machine == "macfrank.rice.edu")*3 + (machine == "star3.local")*2 + (machine == "star4.local")*2 )

InitialDir = /home/sy34/workspace/LowPtMuonAna/bin/
Executable = /home/sy34/workspace/LowPtMuonAna/bin/pairAna.app

GetEnv     = True

Queue 100
