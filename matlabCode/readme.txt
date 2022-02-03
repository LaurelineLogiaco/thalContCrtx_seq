This code was originally created by Laureline Logiaco

Here, the optimization of the cost functions uses CMA-ES (an evolutionary algorithm, http://cma.gforge.inria.fr/)
because in general our cost functions admit several solutions. This is generally less fast (because the optimization is a gradient-free method)
but this method searches for lower minima and gives robust results.

While code for the first figures of the article is already accessible and should be very understandable, we are still in the process of commenting some of the functions, and the full package will be completed soon.

1. fitFnctn_sumOfEigenmodes contains code related to Supp. Fig. S1; which relates to finding a decomposition of any target function into a sum of a few eigenmodes that will constitute a robust solution for the production of the target by the recurrent network. The file to run is MainFitAndPlot.m .

2. 
