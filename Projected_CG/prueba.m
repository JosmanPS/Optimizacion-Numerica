
% Set seed
rng(130056)

n = 10;
cut = 2;
m = round(n^2 / cut);
maxiter = 50;
tol = 1e-6;

G = gallery('poisson',n);
c = rand(n^2, 1);

A = G(1:m, :);
b = rand(m, 1);

x = A \ b;

z = projected_cg(x, G, c, A, b, maxiter, tol)