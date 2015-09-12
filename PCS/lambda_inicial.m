function [lambda] = lambda_inicial(g, A)
    
% -----------------------------------------------------
%
% DESCRIPCIÓN:
% 
% Esta función calcula un valor inicial los parámetros
% de las restricciones planteadas en el Lagrangeano.
% Utilizamos la aproximación de mínimos cuadrados para
% lograr nuestro objetivo.
%
% AUTOR: José Manuel Proudinat Silva
%
% INPUT:
%     - g : El gradiente de la función objetivo
%     - A : El gradiente de las restricciones
% OUTPUT:
%     - lambda : Los multiplicadores de Lagrange
% 
% -----------------------------------------------------

    % Hacemos una descomposición QR para MC
    [Q, R] = qr(A');
    lambda = - Q' * g;
    lambda = R \ lambda;

end
