
    A = rand(100, 10);
    b = rand(100, 1);

    %
    % Método QR
    %
    tic;
    [Q, R] = qr(A);
    x_qr = Q' * b;
    x_qr = R \ x_qr;
    toc_qr = toc;
    clear Q
    clear R

    %
    % Método \
    %
    tic;
    x_back = A \ b;
    toc_back = toc;

    %
    % Método LDL'
    %
    tic;
    K = [  speye(100),    sparse(A)   ;
           sparse(A'),  zeros(10, 10)];
    F = sparse([b; zeros(10,1)]);
    x_ldl = backsolve(K, F);
    x_ldl = x_ldl(101:110);
    toc_ldl = toc;
    