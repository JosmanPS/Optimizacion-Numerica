function [ x,lm,norma,f,i ] = NewtonKKT( nombre, maxiter,tol )
% Joaquin Sanchez 130518 - 0ptimizacion Numerica
%Este metodo resuelve el problema de minimizacion con restricciones
%        min f(x)
%         s.a. c(x) = 0 
% Input: nombre - String del nombre del archivo que se usa .nl en AMPL
%        maxiter - maximo de iteraciones
%        tol - tolerancia

%path(path,'/home/joaxchon/Desktop/Morales'); 

% Concatenamos el nombre con .nl para tener el archivo
nombreAMPL = strcat(nombre, '.nl');

% Usamos AMPL para obtener las funciones gradientes y hessiana
[ x, ~, ~, lm, ~,~ ] = spamfunc(nombreAMPL);

n = length(x);      % numero de variables
m = length(lm);    % numero de restricciones

[ f, c ] = spamfunc ( x, 0 );

[ g, A ] = spamfunc ( x, 1 );

% En el primer paso obtenemos lambda0 por minimos cuadrados
lm = CalculaLambda(g,A);
% gradiente de la Lagrangiana
gL = g - A'*lm;
% Obtenemos la Hessiana de la Lagrangiana
[ W ] = spamfunc ( -lm );
%Un bonito contador
i=0;
% Calculamos la norma
norma=norm(g,2);
normaL=norm(gL,2);
% Loop
while i < maxiter && norma > tol
    % Tenemos la iteracion de la matriz KKT en otra funcion
    [p,lambda]=IteracionNewton(W,A,g,c);
    % actualizamos
    x = x+p;
    lm =  -lambda;
    [ f, c ] = spamfunc ( x, 0 );
    [ g, A ] = spamfunc ( x, 1 );
    [ W ] = spamfunc ( lm );
    norma=norm(g,2);
    gL=g-A'*lm;
    normaL=norm(gL);
    i=i+1;
end
fprintf( ' Nombre del problema                      %s  \n', nombre);
fprintf( ' Numero de variables                    %4i \n', n);
fprintf( ' Numero de restricciones                %4i \n', m);
fprintf( ' Numero maximo de iteraciones            %4i \n', maxiter);
fprintf( ' Tolerancia                               %8.2e \n\n', tol);   
fprintf( ' Objetivo en el ultimo punto             % 21.15e \n', f);
fprintf( ' Norma de las restricciones              % 8.2e \n', norm(c, inf));
fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', normaL );

end

