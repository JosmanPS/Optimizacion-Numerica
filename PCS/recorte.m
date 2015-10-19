function [ alpha, feval ] = recorte(x, f, lm, W, norm_c1, p, mu, clow)
    
    %
    %
    %

    alpha = 1;
    x_alpha = x + alpha * p;    
    [f_alpha, c_alpha] = spamfunc(x_alpha, 0);
    feval = 1;
    c_alpha = c_alpha - clow;
    norm_c1_alpha = norm(c_alpha, 1);
    phi_alpha = f_alpha + mu * norm_c1_alpha;

    phi = f + mu * norm_c1;
    norm_lm = norm(lm, inf);
    D_phi = -p'*W*p - (mu - norm_lm) * norm_c1;
    c1 = 1e-4;
    wolfe = phi + c1 * alpha * D_phi;
    iter = 0;
    
    while phi_alpha > wolfe

        alpha = alpha / 2;
        x_alpha = x + alpha * p;    
        [f_alpha, c_alpha] = spamfunc(x_alpha, 0);
        feval = feval + 1;
        c_alpha = c_alpha - clow;
        norm_c1_alpha = norm(c_alpha, 1);
        phi_alpha = f_alpha + mu * norm_c1_alpha;

        iter = iter + 1;
        
        if iter == 20
            error('Se ha sobrepasado el maximo de iteraciones en el recorte')
        end
     
    end
    
end
