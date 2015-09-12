function [ ] = pcs(nombre, itermax, tol)

  % ----------------------------------------------------------
  %
  %             PROGRAMACIÓN CUADRÁTICA SUCESIVA
  %
  % ----------------------------------------------------------

  % ----------------------------------------------------------
  %
  % DESCRIPCIÓN:
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
  % AUTHOR : José Manuel Proudinat Silva
  %
  % FECHA  : 04/09/15
  %
  % ----------------------------------------------------------

      nombreAMPL = strcat('/home/josmanps/Projects/Optimizacion-Numerica/ampl-models/', nombre, '.nl');

      [ x, xlow, xupp, lm, clow, cupp ] = spamfunc(nombreAMPL);

      n = length(x);       % Número de variables
      m = length(lm);     % Número de restricciones

      % Evaluamos en el punto inicial
      [f, c] = spamfunc(x, 0);
      c = c - clow;
      [g, A] = spamfunc(x, 1);

      % Calculamos los multiplicadores de Lagrange iniciales
      % por medio de Mínimos Cuadrados
      lm = lambda_inicial(g, A);
      
      % Evaluamos el gradiente de la Lagrangeana
      gL = g - A'*lm;
      % Evaluamos la Hessiana de la Lagrangeana en el punto inicial
      [W] = spamfunc(-lm);

      fprintf( ' Nombre del problema                      %s  \n', nombre);
      fprintf( ' Numero de variables                    %4i \n', n);
      fprintf( ' Numero de restricciones                %4i \n', m);
      fprintf( ' Numero maximo de iteraciones            %4i \n', itermax);
      fprintf( ' Tolerancia                               %8.2e \n\n', tol);   
      fprintf( ' Objetivo en el punto inicial            % 21.15e \n', f);
      fprintf( ' Norma de las restricciones              % 8.2e \n', norm(c, inf));
      fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', norm(gL) ...
               );

      %
      % Comenzamos el proceso iterativo de Newton
      %
      iter = 0;
      norm_g = norm(g, 2);
      norm_L = norm(gL, 2);

      while iter < itermax && norm_g > tol
          
          % Calculamos la dirección de Newton para esta iteración
          [p, lm, spd] = Newton(c, g, A, W);

          if not(spd)
              fprintf('\n   *****La inercia es incorrecta***** \n')
              break;
          end
          
          % Actualizamos
          x = x + p;
          lm = -lm;
          [f, c] = spamfunc(x, 0);
          c = c - clow;
          [g, A] = spamfunc(x, 1);
          [W] = spamfunc(lm);
          
          iter = iter + 1;

      end

      %
      % Impresiones finales
      %
      gL = g - A'*lm;
      norm_g = norm(g, 2);
      norm_L = norm(gL, 2);
      
      fprintf( '\n\n\n' )
      fprintf( 'Resultados finales \n' );
      fprintf( '------------------- \n' )
      fprintf( ' Objetivo en el punto final              % 21.15e \n', f);
      fprintf( ' Norma de las restricciones              % 8.2e \n', norm(c, inf));
      fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', norm(gL) ...
               );
      fprintf( ' Número de iteraciones                      %3i \n', iter);

end
  