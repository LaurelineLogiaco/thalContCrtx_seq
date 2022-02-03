This code was originally created by Laureline Logiaco

Here, the optimization of the cost functions uses CMA-ES (an evolutionary algorithm, http://cma.gforge.inria.fr/)
because in general our cost functions admit several solutions. This is generally less fast (because the optimization is a gradient-free method)
but this method searches for lower minima and gives robust results.

While code for the first figures of the article is already accessible and should be very understandable, we are still in the process of commenting some of the functions, and the full package will be completed soon.

1. The folder fitFnctn_sumOfEigenmodes contains code related to Fig. 1 F and Supp. Fig. S1; which relates to finding a decomposition of any target function into a sum of a few eigenmodes that will constitute a robust solution for the production of the target by the recurrent network. The file to run is MainFitAndPlot.m, it illustrates how to use the function cost_func.m (which itself uses EignMd_Approx.m) in order to find complex amplitudes and timescales of a few eigenmodes such that their weighted sum approximates a desired output.

2. The folder half-random_EigenvalueControl contains the function Create IllustratingEignvalueCntrl.m that illustrates how to control eigenvalues with a half random loop, and what the limitations are (related to Fig. 1 E and Supp. Fig. S2).

