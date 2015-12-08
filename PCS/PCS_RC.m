function x = PCS_RC(name, TOL, maxiter)

    nameAMPL = strcat('/home/josmanps/Projects/Optimizacion-Numerica/ampl-models/', name);
    nameAMPL = strcat(nameAMPL, '.nl');

    % NOCEDAL p. 549

    %
    % Parameters
    %
    rho_mu = 0.1;
    etha = 0.1;
    delta = 1;
    % inc = 1.5;
    % dec = 0.8;

    %
    % Compute initial point x
    %
    [ x, ~, ~, ~, clow, ~] = spamfunc(nameAMPL);

    %
    % Compute fk, ck, g, A
    %
    [f, c] = spamfunc(x, 0);
    n_feval = 1;
    c = c - clow;
    [g, A] = spamfunc(x, 1);

    [m, n] = size(A);

    %
    % Compute Lagrange multipliers
    %
    lm = - LDL_ls(A', g);  % TODO : Revisar signo
    mu = norm(lm, inf);

    %
    % Compute gradient and Hessian of Lagrange
    %
    gL = g - A' * lm;
    norm_gL = norm(gL, inf);
    norm_c = norm(c, 1);
    [ W ] = spamfunc(-lm);

    TOL1 = TOL * (1 + norm_gL);
    TOL2 = TOL * (1 + norm_c);

    iter = 0;

    fprintf(['\n iter          f_k             ||c_k||       ||gL_k||        ' ...
            ' mu           delta \n']);
    fprintf(' ---------------------------------------------------------------------------------- ');
    fprintf('\n %3i    %1.11e   %1.5e   %1.5e   %1.5e   %1.5e', ...
            iter, f, norm_c, norm_gL, mu, delta);

    while iter < maxiter

        %
        % Check for optimality
        %

        if norm_gL < TOL1 && norm_c < TOL2
            return;
        end

        %
        % Solve normal subproblem for v and compute r
        %
        [v, ~] = Dogleg(A, c, delta);
        r = A * v + c;

        %
        % Compute p by projected CG
        %
        [p, GC_iter, spd] = PCG_TR(W, g, A, r-c, delta, 1e-6, v);

        %
        % Compute mu
        %
        pWp = p' * W * p;
        if pWp > TOL
            sigma = 0.5;
        else
            sigma = 0;
        end
        %
        norm_cp = norm(c + A*p, 1);
        mu_aux = (g'*p + sigma * pWp) / (norm_c - norm_cp) + 0.1;
        mu = max(mu, mu_aux);

        %
        % Compute rho
        %
        xx = x + p;
        [ff, cc] = spamfunc(xx, 0);
        cc = cc - clow;
        norm_cc = norm(cc, 1);
        [gg, AA] = spamfunc(xx, 1);
        ll = - LDL_ls(AA', gg);    % TODO : Revisar signo
        [WW] = spamfunc(-ll);
        n_feval = n_feval + 1;
        %
        ared = f - ff + mu * (norm_c - norm_cc);
        pred = -p'*g - 0.5 * p'*W*p + mu * (norm_c - norm_cp);
        rho = ared / pred;

        if rho > 0 && ared > 0

            % Update
            x = x + p;
            f = ff;
            c = cc;
            g = gg;
            A = AA;
            W = WW;
            lm = ll;    % TODO : Check correct sign
            gL = g - A' * lm;

            norm_c = norm_cc;
            norm_gL = norm(gL, inf);

            % Choose delta bigger
            delta = 1.5 * delta;

        else

            % Choose smaller delta
            norm_p = norm(p);
            % delta = (1-etha) / (1-rho) * norm_p;
            
            % if delta > 0.5 * norm_p
            %     delta = 0.5 * norm_p;
            % elseif delta < 0.1 * norm_p
            %     delta = 0.1 * norm_p;
            % end
            delta = 0.8 * norm(p);

        end

        iter = iter + 1;

        fprintf('\n %3i    %1.11e   %1.5e   %1.5e   %1.5e   %1.5e', ...
                iter, f, norm_c, norm_gL, mu, delta);

    end

end
                

