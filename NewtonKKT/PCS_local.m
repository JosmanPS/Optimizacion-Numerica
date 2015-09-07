function PCS_local ( nombre, itermax, tol );
%
% directorio en donde reside spamfunc.mexmaci64
%
% path(path,'/Users/jmorales/bin');   
%
% ... Este procedimiento ilustra la comunicacion entre Matlab y AMPL. El
% objetivo es resolver el problema de optimizacion
%
%                    minimizar   f(x)
%                    sujeta a    c(x) = 0.
%  
% utilizando programacion cuadratica sucesiva. AMPL evalua los
% gradientes y Hessianas de f y c por derivacion automatica 
%
%        j-l morales 
%        Departamento de Matematicas
%        ITAM
%        agosto 2015
%
%--------------------------------------------------------------------------
%
nombreAMPL = strcat('/home/josmanps/Projects/Optimizacion-Numerica/ampl-models/', nombre, '.nl');

[ x0, xlow, xupp, lm_0, clow, cupp ] = spamfunc(nombreAMPL);

n = length(x0);      % numero de variables
m = length(lm_0);    % numero de restricciones
% 
% ... evaluar el objetivo y las restricciones en el punto inicial
%
[ f0, c0 ] = spamfunc ( x0, 0 );
%
% ... evaluar el gradiente de f y la Jacobiana de c en el punto inicial
%
[ g0, A0 ] = spamfunc ( x0, 1 );
%
% ... evaluar el gradiente de la Lagrangiana
%
gL = g0 - A0'*lm_0;
%
% ... evaluar la Hessiana de la Lagrangiana en el punto inicial
%
[ W0 ] = spamfunc ( -lm_0 );

fprintf( ' Nombre del problema                      %s  \n', nombre);
fprintf( ' Numero de variables                    %4i \n', n);
fprintf( ' Numero de restricciones                %4i \n', m);
fprintf( ' Numero maximo de iteraciones            %4i \n', itermax);
fprintf( ' Tolerancia                               %8.2e \n\n', tol);   
fprintf( ' Objetivo en el punto inicial            % 21.15e \n', f0);
fprintf( ' Norma de las restricciones              % 8.2e \n', norm(c0, inf));
fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', norm(gL) );
