%
function IPM =  ipm_on ( name, TOL , maxiter );
    
    %--------------------------------------------------------------------------
    %   
    % IPM        An Interior Point Method 
    %            developed by Josman at ITAM.
    %            December, 2015.
    %     
    %            Features:   
    %
    %             * AMPL interface 
    %             * NO fixed variables allowed
    %             * free variables allowed
    %
    %
    % INPUT:
    % ---------
    %   - name    : *.nl file  | AMPL model 
    %   - TOL     : float      | numeric zero tolerance
    %   - maxiter : int        | max number of iterations
    %
    % OUTPUT:
    % ---------
    %   - **
    % 
    %                
    %--------------------------------------------------------------------------


    %--------------------------------------------------------------------------

    %
    % AMPL Values
    %
    point = 'amplpnt';
    fungrad = 'amplstub';
    factor = 0.1;
    maxbar = 10;
    mu0 = 0.1;
    separate = 1;

    %
    % ... open output file
    %
    name = strcat('/home/josmanps/Projects/Optimizacion-Numerica/ampl-models/', name);
    name1 = strcat(name,'.out1'); fout1 = fopen(name1,'w');   
    %
    % ... get initial point from stub; compute initial values: 
    %     obj function and derivatives. 
    %
    stub  = strcat(name, '.nl');

    %
    % INITIAL VALUES
    %
    [ x, s, y ] = feval ( point );
    % 
    n_x = length(x);        % ... number of variables
    n_s = length(s);        % ... number of slacks
    s   = zeros(n_s,1);     %
    %
    [ f, c, g, A, ~ ] = feval (fungrad, x, s, y );
    %
    m   = length(c);        % ... number of c/s (incl. bnds)
    m_e = m - n_s;          % ... number of equality c/s

    thresh  =  1.0d-1;      % ... minimum value for the initial slacks
    %
    s  = Initials ( c, m_e, n_s, thresh );
    e  = ones(n_s, 1);     % ... all ones vector
    %
    % ... set initial barrier, Lagrange multipliers, constraints
    %     compute W with LS Lagrange multipliers
    %
    mu = mu0;

    infeas  = norm(c,2);  
    bar  = f  - mu*log( s )'*e;

    y = Lagrange( A, s, mu, g, e );

    %
    % New initial values
    %
    [ f, c, g, A, W ] = feval (fungrad, x, s, y );
    c_h = c(1:m_e); 
    c_g = c((m_e+1):end);
    A_h = A(1:n_x, 1:m_e);
    A_g = A(1:n_x, (m_e+1):size(A, 2));
    lm_h = y(1:m_e);
    lm_g = y((m_e+1):end);
    inertia_n = n_s + n_x;
    inertia_m = m;
    
    iter = 0;
    in_iter = 0;

    %
    % Logarithmic barrier infeasibility
    %
    infeas  = norm(c,inf);  
    bar  = f  - mu * sum(log( s ));
    gL = g + [A_h A_g] * y;
    norm_gL = norm(gL, inf);
    tol1 = TOL * (1 + norm_gL);
    tol2 = TOL * (1 + infeas);


    %
    % ... print characteristics of the problem 
    %
    Print_head ( fout1, name  , n_x, m  , n_s   , m_e    , factor,...
                     mu0  , maxbar, f  , bar, infeas, separate );

    fprintf(fout1, '\n\n\n   iter        f            f - logbar        || gL ||       || c ||           mu \n');
    fprintf(fout1, '--------------------------------------------------------------------------------------\n');
    fprintf(fout1, '%5i    %1.7e    %1.7e    %1.5e    %1.5e    %1.5e \n', ...
            iter, f, bar, norm_gL, infeas, mu);

    %
    % Start iterative Newton
    %
    while (norm_gL > tol1 || infeas > tol2) && iter < maxiter
        
        %
        % Create the subproblem
        %
        S = spdiags(s, 0, n_s, n_s);
        S_inv = spdiags(1./s, 0, n_s, n_s);

        bar_iter = 0;
        S_obj = -mu * e;
        Ss = S * lm_g + S_obj;
        norm_G = norm([gL; S_inv * Ss], inf);
        obj = -[g; S_obj; c_h; c_g];
        tol3 = mu * max(norm_G, infeas);

        %
        % Solve IPM subproblem
        %
        while (infeas > tol3 || norm_G > tol3) && bar_iter < maxbar

            %
            % Compute the KKT matrix
            %
            L_g = spdiags(lm_g, 0, n_s, n_s);
            WW = [W, sparse(n_x, n_s); sparse(n_s, n_x), S*L_g];
            KKT = [WW, A; A', sparse(m, m)];

            %
            % Solve KKT and check for good inertia
            %
            [d, neg, ran] = backsolve(KKT, obj);
            inertia = (neg==inertia_m && ran==(inertia_n + inertia_m));

            if not(inertia)
                warning('*** Incorrect Inertia ***');
                return;
            end

            %
            % Split and cut the directions
            %
            dx = d(1:n_x);
            ds = S * d((n_x + 1):(n_x + n_s));
            dy = d((n_x + n_s + 1):end) - y;
            dlg = d((n_x + n_s + m_e + 1):end) - lm_g;

            alpha_s = step(s, ds, 0.995);
            alpha_lg = step(lm_g, dlg, 0.995);

            %
            % Update variables
            %
            x = x + alpha_s * dx;
            s = s + alpha_s * ds;
            y = y + alpha_lg * dy;

            %
            % Compute and update the subproblem
            %
            S = spdiags(s, 0, n_s, n_s);
            S_inv = spdiags(1./s, 0, n_s, n_s);

            [ f, c, g, A, W ] = feval (fungrad, x, s, y );
            c_h = c(1:m_e); 
            c_g = c((m_e+1):end);
            A_h = A(1:n_x, 1:m_e);
            A_g = A(1:n_x, (m_e+1):size(A, 2));
            lm_h = y(1:m_e);
            lm_g = y((m_e+1):end);

            gL = g + [A_h A_g] * y;

            S_obj = -mu * e;
            Ss = S * lm_g + S_obj;

            obj = -[g; S_obj; c_h; c_g]; % Lado derecho del sistema KKT.
            infeas = norm(c, inf);
            bar = f - mu * sum(log( s ));
            norm_gL = norm(gL, inf);
            norm_G = norm([gL; S_inv * Ss], inf);

            bar_iter = bar_iter + 1;
            in_iter = in_iter + 1;

            fprintf(fout1, '%5i    %1.7e    %1.7e    %1.5e    %1.5e    %1.5e \n', ...
                    in_iter, f, bar, norm_gL, infeas, mu);

        end

        mu = mu * factor;
        iter = iter + 1;

    end

end
%
% -------------------------------------------------------------------------
%


function s = Initials ( c, n_e, n_s, thresh );
    zero = 0.0d0;
    %
    % ... This subroutine computes initial values for the slacks 
    %     associated with inequality constraints c(x) <= 0 
    %
    %     if  c(x0) <=0   AND   c(x0) < thresh   
    %         s0 =  c(x0)
    %     else 
    %         s0 =  thresh
    %     end 
    %
    s = zeros(n_s,1);
    for i=1:n_s
        j  = n_e + i;
        cj = c(j);
        if cj <= zero & cj < -thresh
            s(i) = -cj;
        else
            s(i) = thresh;
        end
    end

end
%
% -------------------------------------------------------------------------
%


function print = Print_head ( fout , name  , n_x   , m       , n_s, ...
                              m_e  , factor, mu0   , maxbar, ...
                              f    , bar   , infeas, separate );
    %
    % ... this routine prints initial values                          
    %
    if separate == 0
       SEP = 'Equal   ';
    else
       SEP = 'Separate';
    end
    %
    fprintf(fout,'\n An Interior Point Method solving problem: %s \n', name);
    fprintf(fout,' Problem dimensions: n_x = %4i;\tm   = %4i            \n', n_x, m);
    fprintf(fout,'                     n_s = %4i;\tm_e = %4i            \n', n_s, m_e);
    fprintf(fout,'                                                      \n');
    fprintf(fout,' Parameters:         factor ................. %7.1e   \n', factor);
    fprintf(fout,'                     mu0 .................... %7.1e   \n', mu0);
    fprintf(fout,'                     maxbar ................. %4i     \n', maxbar);
    fprintf(fout,'                     Step lengths ........... %s      \n', SEP );
    fprintf(fout,'                                                     \n');
    fprintf(fout,' Objective   ............... % 14.8e   \n', f);
    fprintf(fout,' Barrier     ............... % 14.8e   \n', bar);
    fprintf(fout,' ||c||_2     ............... % 14.8e   \n', infeas);

end
%
%--------------------------------------------------------------------------
%


function y = Lagrange(A, s, mu, g, e)

    % Least squares Lagrangian multipliers

    [n, m] = size(A);
    
    OBJ = A' * [-g; mu * e];
    OBJ = [ sparse(n, 1)   ;
               OBJ        ];
    
    K = [  speye(n, n) ,       A        ;
              A'       ,  sparse(m, m) ];
    y = backsolve(K, OBJ);
    y = y((n + 1):end);
    n_s = length(s);
    n_h = m - n_s;

    if n_s ~= 0

        y_y = y((n_h + 1):(n_h + n_s));
        ind = (y_y <= 0);
        aux = mu ./ s;
        y_y(ind) = min(1e-3, aux(ind));

        y((n_h + 1):(n_h + n_s)) = y_y;

    end

end
%
%--------------------------------------------------------------------------
%
