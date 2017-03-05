# Wave-Equation-Parallelization. Meant to be run on the SciNet GPC.

Here is code that parallelizes an implementation of a wave equation evolution to arbitrary nodes. 

Visualization of wave equation evolution on 1 node: https://www.youtube.com/watch?v=-qsKpPFmc3E

Visualisation of parallelizing the wave equation implementation with the 3 nodes: https://www.youtube.com/watch?v=oXrn8Y5xI4k

There are folders for the plots, the report which provides analysis, the git log, and 2 folders with the main code (wave1d.cc). The code in "with-IO" has the IO ready so that you can test the visuals. The code in "no-IO" is the code I used for the scaling as it has the IO commented out. 

Here is a scaling analysis, you can see more figures in the plots folder.
![alt tag](http://i.imgur.com/1FVDwz8.png)
