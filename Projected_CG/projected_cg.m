function [ x, iter, feval ] = projected_cg(x, G, c, A, b, maxiter, tol)
    
    % ----------------------------------------------------------
    %
    %                GRADIENTE CONJUGADO PROYECTADO
    %
    % ----------------------------------------------------------

    % ----------------------------------------------------------
    %
    % DESCRIPCIÓN:
    %
    % El objetivo es resolver el problema:
    %
    %                    min    1/2 p' G p + p' g
    %                    s.a.   A p + c = 0       
    % 
    % Donde A es de m x n, con m < n.
    % Para resolver este problema se utiliza la teoría de Gradiente
    % Conjugado Proyectado.
    % Como referencia de la teoría detrás del método se puede tomar:
    %     J.Nocedal & S. Wright, "Numerical Optimization", cp. 16,
    %     2nd Edition, Springer, 2006
    %
    % AUTOR : José Manuel Proudinat Silva
    %
    % COMMITS :
    %    
    %     * 2015-10-23 | Método inicial
    %           Se hizo un método inicial que puede ir mejorando mientras
    %           se explore más teoría.
    %
    % ----------------------------------------------------------

    % Calculamos x un punto inicial factible
    % x = A \ b;
    [m, n] = size(A);
    iter = 0;
        
    % Calculamos los residuos y direcciones proyectadas
    r = G * x + c;
    % Generamos la matriz que nos ayudará a encontrar la solución proyectada
    K = [ speye(n),       A'      ;
          A       ,  sparse(m, m)];
    k = [r; sparse(m, 1)];
    g = backsolve(K, k);
    g = g(1:n);
    d = -g;
    norm_g = norm(g, inf);
    
    while iter <= maxiter && norm_g > tol;

        Gd = G * d;
        dGd = d' * Gd;
        rg = r' * g;
        alpha = rg / dGd;

        x = x + alpha * d;
        r = r + alpha * Gd;

        k = [r; sparse(m, 1)];
        g = backsolve(K, k);
        g = g(1:n);
        d = -g;
        norm_g = norm(g, inf);

        beta = (r' * g) / rg;
        d = -g + beta * d;

    end
    
end

