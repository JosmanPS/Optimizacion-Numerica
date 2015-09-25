function x = pre_gc( H, g)

% =========================================================
%
% Este programa encuentra la solucion (x) del problema:
%				Hx = g
% a traves del metodo de gradiente conjugado precondicionado.
% Usamos la versión más simple de la factorización de 
% Cholesky incompleta, en este caso.
%
% 24 Septiembre 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% Entrada:
% 	- H : matriz del problema (Hessiana)
% 	- g : vector de respuestas (gradiente)
%
% Salida:
% 	- x : el vector solucion al problema 
%
% =========================================================

    % Calculamos los valores iniciales
    n = length(g);
    maxit = 2*n;
    x = zeros(n,1);
    r = H * x - g;
    iter = 0;

    % Precondicionamos
    opts.michol = 'on';
    L = ichol(H, opts);
    y = L \ r;
    y = L' \ y;
    d = -y;
    rTy = r' * y;
    
    % Guardamos los resultados para reducir operaciones
    rTr = r' * r;
    Hd = H * d;
    dHd = d' * Hd;
    TOL_r0 = 1e-16 * rTr;
    TOL2 = 1e-6;

    while(iter < maxit && rTr > TOL_r0 && dHd > TOL2)

        % Actualizamos los valores
        alpha = rTy / dHd;
        x = x + alpha * d;
        r = r + alpha * Hd;
        y = L \ r;
        y = L' \ r;
        Beta = (r' * y) / rTy;
        d = -y + Beta * d;

        % Guardamos los resultados para reducir operaciones
        rTr = r' * r;
        rTy = r' * y;
        Hd = H * d;
        dHd = d' * Hd;
        
        % Siguiente iteracion
        iter = iter + 1;

    end

end
