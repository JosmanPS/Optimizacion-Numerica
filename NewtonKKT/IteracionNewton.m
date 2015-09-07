function [ p,lambdaop ] = IteracionNewton( W,A,g,c )

% Este programa calcula el avance en x (p) y la nueva -lambda
% Considera la factorizacion LDL de la matriz de Karush-Kuntucker
% Input: W- La Hessiana de la Lagrangiana
%        A = la matriz de gradientes de las restricciones
%        g  = el gradiente de la funcion f
%        c = la evaluacion de las restricciones 
b = -[g; c];
[w1, w2] = size(W);
[m, n] = size(A);
K = [W, A'; A, zeros(m)];
[L, D, P] = ldl(K);
%Tengo que checar si tengo que usar la P (yo creo que si)
z = P * L \ b;
x = (D * L' * P')\z;
p = x(1:w2);
lambdaop = x(w2+1:w2+m);




end

