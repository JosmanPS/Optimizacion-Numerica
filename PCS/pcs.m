function [ x, n, m, f, iter, feval, time, spd ] = pcs(nombre, itermax, tol, verbose)

  % ----------------------------------------------------------
  %
  %             PROGRAMACIÓN CUADRÁTICA SUCESIVA
  %
  % ----------------------------------------------------------

  % ----------------------------------------------------------
  %
  % DESCRIPCIÓN:
  %
  % Este procedimiento ilustra la comunicacion entre Matlab y AMPL.
  % El objetivo es resolver el problema de optimizacion:
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
      tic;
      n = length(x);       % Número de variables
      m = length(lm);     % Número de restricciones

      % Evaluamos en el punto inicial
      [f, c] = spamfunc(x, 0);
      feval = 1;
      c = c - clow;
      norm_c1 = norm(c, 1);
      norm_c_inf = norm(c, inf);
      [g, A] = spamfunc(x, 1);

      % Calculamos los multiplicadores de Lagrange iniciales
      % por medio de Mínimos Cuadrados
      % lm = lambda_inicial(g, A);
      lm = -A' \ g; 

      % Agregamos los valores iniciales de la funci'on de m'erito
      rho = 1e-1;
      mu = norm(lm, inf);
      
      % Evaluamos el gradiente de la Lagrangeana
      gL = g - A'*lm;
      norm_gL = norm(gL, inf);
      % Evaluamos la Hessiana de la Lagrangeana en el punto inicial
      [W] = spamfunc(lm);

      if verbose
          fprintf( ' Nombre del problema                      %s  \n', nombre);
          fprintf( ' Numero de variables                    %4i \n', n);
          fprintf( ' Numero de restricciones                %4i \n', m);
          fprintf( ' Numero maximo de iteraciones            %4i \n', itermax);
          fprintf( ' Tolerancia                               %8.2e \n\n', tol);   
          fprintf( ' Objetivo en el punto inicial            % 21.15e \n', f);
          fprintf( ' Norma de las restricciones              % 8.2e \n', norm_c_inf);
          fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', norm_gL);
      end      
          
      %
      % Comenzamos el proceso iterativo de Newton
      %
      iter = 0;

      if verbose
          fprintf('\n ************************************************** \n');
          fprintf(['\n iter          f_k             ||c_k||       ||gL_k||        ' ...
                   ' mu           alpha \n']);
          fprintf(' ---------------------------------------------------------------------------------- ');
          fprintf('\n %3i    %1.11e   %1.5e   %1.5e   %1.5e   %1.5e', ...
                  iter, f, norm_c_inf, norm_gL, mu);
      end
          
      % Calculamos la toleracncia relativa
      tol_gL = tol * (1 + norm_gL);
      tol_c = tol * (1 + norm_c_inf);
      
      while iter < itermax && (norm_gL > tol_gL || norm_c_inf > tol_c)
          
          % Calculamos la dirección de Newton para esta iteración
          [p, dlm, spd] = Newton(c, g, A, W);

          if not(spd)
              fprintf('\n\n   *****La inercia es incorrecta***** \n')
              break;
          end
          
          % Calculamos los pasos
          dlm = -dlm;
          dlm = dlm - lm;
          [alpha, rec_feval] = recorte(x, f, lm, W, norm_c1, p, mu, clow);
          feval = feval + rec_feval;
          
          % Actualizamos
          x = x + alpha * p;
          lm = lm + alpha * dlm;
          [f, c] = spamfunc(x, 0);
          feval = feval + 1;
          c = c - clow;
          [g, A] = spamfunc(x, 1);
          gL = g - A'*lm;
          norm_gL = norm(gL, inf);
          [W] = spamfunc(-lm);

          % Calculamos mu_k
          pWp = p' * W * p;
          sigma = 0;
          if pWp > 0
              sigma = 0.5;
          end
          gp = g' * p;
          norm_c1 = norm(c, 1);
          norm_c_inf = norm(c, inf);
          mu_aux = gp + sigma * pWp;
          mu_aux = mu_aux / ( (1-rho) * norm_c1 );
          mu = max(mu, mu_aux);
          
          iter = iter + 1;

          if verbose
              fprintf('\n %3i    %1.11e   %1.5e   %1.5e   %1.5e   %1.5e', ...
                      iter, f, norm_c_inf, norm_gL, mu, alpha);
          end
          
      end
      fprintf('\n\n')
      time = toc;
      %
      % Impresiones finales
      %
      gL = g - A'*lm;
      norm_g = norm(g, 2);
      norm_L = norm(gL, 2);

      if verbose
          fprintf( '\n\n\n' )
          fprintf( 'Resultados finales \n' );
          fprintf( '------------------- \n' )
          fprintf( ' Objetivo en el punto final              % 21.15e \n', f);
          fprintf( ' Norma de las restricciones              % 8.2e \n', norm(c, inf));
          fprintf( ' Norma del gradiente de la Lagrangiana    %8.2e \n', norm_gL ...
                   );
          fprintf( ' Número de iteraciones                      %3i \n', iter);
      end
          
end
  