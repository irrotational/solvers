# solvers
Simple python examples of the Newton-Raphson and Halley methods, with test functions. To run, simply type:

	python3 solvers.py

To try your own test function, replace 'mytestfunc' in the python script with your own function, which should accept a single argument 'x' and return a single value (the value of the function evaluated at 'x').

To change the range over which roots are searched, change the 'min' and 'max' variables as appropriate. Note that due to the nature of the algorithms (that is; they start at a seed point and then 'move around' to find the root), it is perfectly possible for a root to be found outside of the specified [min,max] range.

To use the Halley method istead of the default Newton-Raphson, pass the 'solver' argument to the 'find_multiple_roots' function as follows:

	find_multiple_roots(mytestfunc,min,max,solver='halley')
