# solvers
Simple python examples of the Newton-Raphson and Halley methods, with test functions. To run, simply type:

	python3 solvers.py

To try your own test function, replace 'mytestfunc' in the python script with your own function, which should accept a single argument 'x' and return a single value (the value of the function evaluated at 'x').

To change the range over which roots are searched, change the 'min' and 'max' variables as appropriate.

To use the Halley method istead of the default Newton-Raphson, pass the 'solver' argument to the 'find_multiple_roots' function as follows:

	find_multiple_roots(mytestfunc,min,max,solver='halley')
