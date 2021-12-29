import numpy as np
import matplotlib.pyplot as plt
from scipy.special import airy,jv,spherical_jn

##################################################################################################
# Derivative functions

def first_derivative(f,x,step_size=1e-6):
	""" Computes the first-order derivative of the callable 'f'
	    using the finite-differences step size 'step_size'. """

	deriv = ( f(x+step_size) - f(x) ) / step_size

	return deriv

def second_derivative(f,x,step_size=1e-6):
	""" Computes the second-order derivative of the callable 'f'
	    using the finite-differences step size 'step_size'. """

	deriv = ( first_derivative(f,x+step_size) - first_derivative(f,x) ) / step_size

	return deriv

##################################################################################################
# Solvers

def newton_raphson(f,initial_guess,prec=1e-6,max_steps=1000000):
	""" Finds the first root of the supplied function 'f'
	    (which must be a callable that returns a single value),
	    using the Newton-Raphson method, given the starting initial
	    guess 'initial_guess', to the precision 'prec'. """

	x = initial_guess
	precision_satisfied = False
	step = 0
	while (not precision_satisfied and step<max_steps):
		if (step == max_steps-1):
			break
		fnc = f(x)
		f_p = first_derivative(f,x)
		x_new = x - fnc/f_p
		delta = x_new - x
		if ( abs(delta) < prec ):
			precision_satisfied = True
		x = x_new
		step += 1
	
	if ( step == (max_steps-1) ):
		print("WARNING: Did not find a root in %d steps. Returned value of x is not a root." % (max_steps))
		return x,False
	else:
		return x,True

def halley(f,initial_guess,prec=1e-6,max_steps=1000000):
	""" Finds the first root of the supplied function 'f'
	    (which must be a callable that returns a single value),
	    using Halley's method, given the starting initial guess
	    'initial_guess', to the precision 'prec'. """

	x = initial_guess
	precision_satisfied = False
	step = 0
	while (not precision_satisfied and step<max_steps):
		fnc = f(x)
		f_p = first_derivative(f,x)
		f_pp = second_derivative(f,x)
		x_new = x - ( 2*fnc*f_p ) / ( 2*f_p**2 - fnc*f_pp )
		delta = x_new - x
		if ( abs(delta) < prec ):
			precision_satisfied = True
		x = x_new
		step += 1

	if ( step == (max_steps-1) ):
		print("WARNING: Did not find a root in %d steps. Returned value of x is not a root." % (max_steps))
		return x,False
	else:
		print("Took %d steps" % (step))
		return x,True

##################################################################################################

def find_multiple_roots(f,min,max,max_guesses=100,solver='newton-raphson'):
	""" Bootstraps the supplied solving function 'solver' using a
	    random guess each time, in order to potentially find
	    multiple roots of the supplied function 'f'. The guesses
	    all lie in the range ['min','max'], and this function will
	    iterate until 'max_guesses' guesses have been exhausted.

	    There is no guarantee that all roots will be found, especially
	    if 'max_guesses' is small. """

	guess_count = 0
	roots = []
	while (guess_count < max_guesses):
		guess = np.random.random() * (max-min) + min
		if (solver == 'newton-raphson'):
			x_root,success_flag = newton_raphson(mytestfunc,guess)
		elif (solver == 'halley'):
			x_root,success_flag = halley(mytestfunc,guess)
		else:
			print('ERROR: Solver not recognised!')
			exit()

		# Check whether this is truly a new root, to some tolerance
		new_root = True
		if ( len(roots) == 0 ): # If we've not yet found any roots, it must be a new root
			new_root = True
		else: # Otherwise, check that this new root differs from all others by some tolerance
			for root in roots:
				diff = root - x_root
				if ( abs(diff) < 1e-3 ):
					new_root = False
		if ( success_flag and new_root ):
			roots.append(x_root)

		guess_count += 1
	
	num_roots_found = len(roots)

	print("Found %d roots in %d guesses" % (num_roots_found,max_guesses))

	return roots
	
##################################################################################################
# Test functions

def mytestfunc(x):
	""" Some weird complicated test function I came up with.
		I count 10 roots in the range [-10,10]. """
		
	func = x**2 * (x**2 * jv(3,x) - 1) * airy(x)[0] + 0.01 * ( np.cosh(x) - 0.1*np.exp(x) )
	return func

##################################################################################################
# Run solver

# NOTE: The solvers may find roots outside of the range [min,max] due to the nature of the algorithms

min,max = -10,10
grid = np.linspace(min,max,1000)
func_grid = mytestfunc(grid)

roots = find_multiple_roots(mytestfunc,min,max)
print(roots)
zeroes = [ 0 for x in roots ]

plt.plot(grid,func_grid)
plt.scatter(roots,zeroes)
plt.show()
