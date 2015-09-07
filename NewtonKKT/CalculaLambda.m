function [ lambda ] = CalculaLambda( g,A )

% Este programa calcula los multiplicadores de Lagrange iniciales 
% usando minimos cuadrados por la factorizacion Qr

    A = A';
    [Q, R] = qr(A);
    z = (R*R') \ g;
    lambda = R' * Q * z;

end

