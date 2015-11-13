c = [-2; -0.5; 0; 0; 0];
Q = sparse(5, 5);
A = [ 1   1   0  1  0 ;
      2   0   0  0  1 ;
      1  8/3  1  0  0 ];
b = [2; 3; 4];
tol = 1e-8;
maxiter = 100;

IPM_QP_plot(Q, c, A, b, tol, inf)