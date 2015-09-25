n = 100;
name = 'poisson';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;


fprintf( 'Problema & n & GC & GCP & pgc & chol \\\\ \\hline \n');
fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 300;
name = 'poisson';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 500;
name = 'poisson';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1000;
name = 'lehmer';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
A = sparse(A);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 3000;
name = 'lehmer';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
A = sparse(A);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 5000;
name = 'lehmer';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
A = sparse(A);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 500;
name = 'toeppd';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
A = sparse(A);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1000;
name = 'toeppd';

% Sistema cuya solucion es el vector de unos
A = gallery(name,n);
A = sparse(A);
n = length(A);
b = A*(ones(n,1));

tic;
x1 = grad_conj(A,b);
t1 = toc;

tic;
x2 = pre_gc(A,b);
t2 = toc;

tic;
[x3, ~] = pcg(A,b);
t3 = toc;

tic;
L = chol(A,'lower');
x4 = L \ b;
x4 = L' \ x4;
t4 = toc;

fprintf( '%s & %4i & %f & %f & %f & %f \\\\ \n', name, n, t1, t2, t3, t4);
