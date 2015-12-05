function [d, neg, ran] = backsolve(A, b)

    [L, D, P, S, neg, ran] = ldl(A);
    
    %
    % Solve the problem
    %
    d = P'*(S*b);
    d = L\d;
    d = D\d;
    d = L'\d;
    d = P'\d;
    d = S*d;
    