function d = backsolve(A, b)

    [L,D,P,S] = ldl(A);
    
    d = P'*(S*b);
    d = L\d;
    d = D\d;
    d = L'\d;
    d = P'\d;
    d = S*d;