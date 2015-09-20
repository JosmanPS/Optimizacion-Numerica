function [p, lm, spd] = Newton(c, g, A, W)
    
% -----------------------------------------------------
%
% DESCRIPCIÓN:
% 
% Esta función resuelve el problema de Newton para
% calcular la dirección de descenso para las variables
% de decisión y los multiplicadores de Lagrange.
%
% AUTOR: José Manuel Proudinat Silva
%
% INPUT:
%     - c  : Las restricciones
%     - g  : El gradiente de la función objetivo
%     - A  : El gradiente de las restricciones
%     - W  : La Hessiana de la Lagrangeana
%
% OUTPUT:
%     - p  : La dirección de descenso
%     - lm : Los multiplicadores de Lagrange
% 
% -----------------------------------------------------

    % Construimos el sistema de KKT
    [m, n] = size(A);
    K = [ W, A'; A, zeros(m)];
    K = sparse(K);
    b = -[g; c];

    % Resolvemos el sistema de KKT
    [L, D, P, S, neg, ran] = ldl(K);
    b = P' * S * b;
    z = L \ b;
    z = D \ z;
    z = L' \ z;
    z = S * P * z;

    % Separamos las direcciones de descenso
    [~, w] = size(W);
    p = z(1:w);
    lm = z((w+1):(w+m));

    % Verificamos la inercia de la matriz
    spd = (neg==m && ran==n+m);

end
