function [ x, iter, norm_rg ] = PCG_TR(G, c, A, b, delta, tol, x)
    
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
    %                    min    1/2 p' G p + c' g
    %                    s.a.   A p = b    
    %                           ||p|| <= delta   
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
    % ----------------------------------------------------------

    %
    % Punto inicial
    %
    [m, n] = size(A);
    K = [ speye(n),       A'      ;
          A       ,  sparse(m, m)];

    if nargin < 7
        k = [sparse(n, 1); b];
        x = backsolve(K, k);
        x = x(1:n);
    end
        
    % Calculamos los residuos y direcciones proyectadas
    r = G * x + c;
    k = [r; sparse(m, 1)];
    g = backsolve(K, k);
    g = g(1:n);
    d = -g;

    norm_g = norm(g, inf);
    rg = r' * g;
    norm_rg = abs(rg);
    % sqrt_r0_y0 = sqrt(rg);
    % err_rel = 1;
    tol = tol * abs(rg);
    
    maxiter = 2000;
    iter = 0;
    
    while iter <= maxiter &&  norm_rg > tol;

        Gd = G * d;
        dGd = d' * Gd;
        alpha = rg / dGd;
        
        %
        % Caso de inercia incorrecta
        %

        % TODO : Verficar caso
        if dGd < 0
            break;
        end

        xx = x + alpha * d;

        % TODO : Verificar caso
        if norm(xx) >= delta
            a = d' * d * alpha^2;
            b = 2 * alpha * d' * x;
            c = x' * x - delta^2;

            tau = (-b + sqrt(b^2 - 4 * a * c)) / (2 * a);
            if imag(tau) ~= 0
                tau = 0;
            end

            tau_aux = (-b - sqrt(b^2 - 4 * a * c)) / (2 * a);
            if imag(tau_aux) ~= 0
                tau_aux = 0;
            end

            if abs(tau - 0.5) > abs(tau_aux - 0.5)
                tau = tau_aux;
            end

            tau = alpha * tau;
            x = x + tau * d;
            return;
        end

        x = x + alpha * d;
        r = r + alpha * Gd;

        k = [r; sparse(m, 1)];
        g = backsolve(K, k);
        g = g(1:n);
        norm_g = norm(g, inf);

        beta = (r' * g) / rg;
        d = -g + beta * d;
        rg = r' * g;
        norm_rg = abs(rg);
        % err_rel = sqrt(rg) / sqrt_r0_y0;

        iter = iter + 1;

    end
    
end
