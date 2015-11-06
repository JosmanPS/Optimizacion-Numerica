
% Set seed
rng(130056)

n = 100;
m = 100;
maxiter = 400;
tol = 1e-6;

G = gallery('poisson', n);
A = G(1:m, :);
n = n^2;

matriz = [ G,     A'       ;
           A, sparse(m, m)];
sol = matriz * ones(n + m, 1);
c = -sol(1:n);
b = sol((n+1):(n+m));

[Q, R] = qr(A');
x = R' \ b;
x = Q * x;

tic;
[z, iter] = projected_cg(x, G, c, A, b, maxiter, tol);
toc