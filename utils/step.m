function alpha = step(x, dx, tau)

% ------------------------------------------------------
%
% Esta función calcula el paso 'alpha', tal que:
%       x = x + alpha * dx 
% se mantenga positivo.
%
% ------------------------------------------------------
    
    neg = (dx < 0.0d0);
    
    if not(sum(neg))
        alpha = tau;
    else
        x = x(neg);
        dx = dx(neg);
        dx = x ./ dx;
        alpha = tau * min( -dx );
        alpha = min(1.0d0, alpha);
    end
        
end
    