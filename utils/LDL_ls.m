function x = LDL_ls(A, b)

    %
    % Symmetric method for Least Squares
    %
    
    [m, n] = size(A);

    K = [  speye(m)  ,       A        ;
              A'     ,  sparse(n, n) ];
    k = [sparse(m, 1); A' * b];

    x = backsolve(K, k);
    x = x((m + 1):end);

end