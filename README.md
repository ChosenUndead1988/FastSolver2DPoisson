# FastSolver2DPoisson
Geometric multigrid solver for the homogeneous Poisson equation on a rectangle using a fourth order compact scheme. 

Consider the Poisson equation 
u_{xx} + u_{yy} = f, (x, y) \in D = [a_x, b_x] \times [a_y, b_y]
where 
u(x, b_x) = 0 = u(x, b_y) when x \in [a_x, b_x]
u(a_x, y) = 0 = u(a_y, y) when y \in [a_y, b_y]. 

Consider the uniform step sizes 
dx = (b_x - a_x)/(N_x + 1) 
dy = (b_y - a_y)/(N_y + 1)

and the coordinates 
x_i = a_x + i*dx, for i = 0, 1, 2, ..., N_x + 1
y_j = a_y + j*dy, for j = 0, 1, 2, ..., N_y + 1

where u(x_i, y_j) = u_{i,j}. The fourth order compact scheme (FOCS) satisfies 

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

