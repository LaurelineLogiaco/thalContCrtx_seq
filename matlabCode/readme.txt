This code was originally created by Laureline Logiaco

Here, the optimization of the cost functions uses CMA-ES (an evolutionary algorithm, http://cma.gforge.inria.fr/)
because in general our cost functions admit several solutions. This is generally less fast (because the optimization is a gradient-free method)
but this method searches for lower minima and gives robust results.

While some is already accessible and should be very understandable, we are still in the process of commenting some of the functions, and the full package will be completed by the end of June. 
