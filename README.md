CASTEP Profiling Visualisation Tool
-----------------------------------

This tool attempts to aid in the quick visualisation of Castep profiling files by graphing the time spent in various routines.

Example usage
------------
Plots the first example profile 
      ./pgraph.py example-profiles/qqq.0001.profile

Plots more data from the first example profile with a total module time cutoff of 5% of the total
      ./pgraph.py example-profiles/qqq.0001.profile -c 5