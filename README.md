# FastSolver2DPoisson
Geometric multigrid solver for the homogeneous Poisson equation on a rectangle using a fourth order compact scheme.

%%%%%%%%%%%%%%%%%%%%%%%% Shorth Description of fourth order compact scheme %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Consider the Poisson equation 
u_{xx} + u_{yy} = f, (x, y) \in D = [a_x, b_x] \times [a_y, b_y]
where 
u(x, b_x) = 0 = u(x, b_y) when x \in [a_x, b_x]
u(a_x, y) = 0 = u(a_y, y) when y \in [a_y, b_y]. 

Consider the uniform step sizes 
dx = (b_x - a_x)/(N_x + 1) 
dy = (b_y - a_y)/(N_y + 1)

and the coordinates 
x_i = a_x + i*dx, for i = 0, 1, 2, ..., N + 1
y_j = a_y + j*dy, for j = 0, 1, 2, ..., N + 1

where u(x_i, y_j) = u_{i,j}. 

The fourth order compact scheme (FOCS) satisfies 

m1*( u_{i + 1, j + 1} + u_{i + 1, j - 1} + u_{i - 1, j + 1} + u_{i - 1, j - 1}) 
  + m2*(u_{i, j + 1} + u_{i, j - 1}) 
  + m3*( u_{i + 1, j} + u_{i - 1, i}) - m4*u_{i,j} 
  = (dx*dx/2.0)*( 8.0*f_{i,j} + f_{i + 1, j} + f_{i - 1, j} + f_{i, j + 1} + f_{i, j - 1}) 

where the coeffcients 

m1 = (1 + \lambda^2)/2.0 
m2 = 5*\lambda^2 - 1
m3 = 5 - \lambda^2 
m4 = 10*(1 + \lambda^2)
 
and \lambda = dx/dy

This results in a sparse linear system L_h u = R_h f = f_R

%%%%%%%%%%%%%%%%%%%%%%%% HOW to set the multigrid method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Our multigrid settings are stored in FOCS2DSettings.f90 
We will enumerate the parameters defined in in this module  

MULTICOLOR_SMOOTHER is a logical variable which defines which version 
of Gauss Seidels is utilized. If its set to true, use the four-color 
version of Gauss-Seidel. Otherwise use (pointwise) lexicographic 
Gauss-Seidel. 

NUMBER_AUXILIARY_GRIDS and LEADING_FACTOR controls the number of grids 
points N. Since geomptric multigrid is a "divide and conquer" algorithm 
we choose  
N = LEADING_FACTOR*2^(NUMBER_AUXILIARY_GRIDS) - 1

ORDER_RESTRICTION_OPERATOR is the order of the restriction operator 
which transfers a grid functions from the fine grid to a coarse grid. 
If set to 2 use Full Weighting operator. If set to 0 use injection. 
Full weighting is treated as the default value.

ORDER_INTERPOLATION_NESTED_ITERATION is the order of the 
prolongation operator. If set to 
- 1 use bi-linear interpolation
- 2 use bi-quadratic interpolation 
- 3 use bi-cubic interpolation 
- 4 use bi-quartic interpolation 

NUM_SWEEPS_DESCENT_MULTIGRID is the number applications of the 
Gauss-Seidel along the "descending portion" of the multigrid 
algorithm. For efficiency set the number ie. 1 or 2. 

NUM_SWEEPS_ASCENT_MULTIGRID is the number of applications of the 
Gauss-Seidel along the "ascending portion" of the multigrid 
algorithm. For efficiency set the number ie. 1 or 2.

NUM_CYCLES_NESTED_ITERATION is the number is iteration per cycle 
for the nested iteration. Nested iteration is a procedure for providing 
a very good initial guess to iterative solver by solving the Poisson 
equation on coarse grids. 

MAXIMUM_CYCLES_MULTIGRID sets a hard cap on the number iterations of multigrid 

CYCLE_TYPE_MULTIGRID chooses the cycle type of the multigrid method. 
Can be either a 'V' cycle or an 'F' cycle. 

ABSOLUTE_TOLERANCE_MULTIGRID and RELATIVE_TOLERANCE_MULTIGRID determines the 
termination criterion of the multigrid iteration. If the residual 
of the i^th iteration is r_h^k = f_R - L_h u^k, then 
the termination critria is 
| r_h^k | \leq RELATIVE_TOLERANCE_MULTIGRID*|r_h^0| + ABSOLUTE_TOLERANCE_MULTIGRID

%%%%%%%%%%%%%%%%%%%%%%%% Test Problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A sample test problem is listed in TestProblems/TestProblem1.f90.
Its given domain is D = [0, 2] \times [0,1] (defined in the module) 
with the test solution u(x,y) = sin(10.0*pi*x)*sin(2.0*pi*y). 

in the Directory containing the makefile makeTest1FOCS.mk run the following commands 
make -f makeTest1FOCS.mk 
./exec 

To delete all the executable files 
make -f makeTest1FOCS.mk clean



