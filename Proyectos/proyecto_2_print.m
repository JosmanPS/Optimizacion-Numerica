function x = proyecto_2_print(ampl_model)

    % ----------------------------------------------------------------------
    %
    % Esta función imprime los resultados del problema evaluado con
    % Programación Cuadrática Sucesiva.
    %
    % Josman, 2015
    %
    % ----------------------------------------------------------------------
    
    [ x, n, m, f, iter, feval, time, spd ] = pcs(ampl_model, 50, 1e-6, 0);

    if spd
        spd = '\cmark';
    else
        spd = '\xmark';
    end
    
    fprintf(['%s  &  %5i  &  %5i  &  %1.14e  &  %3i  &  %3i  &  %1.2e  %s \\\\ ' ...
             ''], ampl_model, n, m, f, iter, feval, time, spd);