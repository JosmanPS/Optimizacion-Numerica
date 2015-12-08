function [ p, type ] = Dogleg ( A, c, Delta );
% ... 
% ... This procedure computes an approximate solution of the problem
%
%            minimize    || Ap + c ||^2
%           subject to   ||p|| <= 0.8Delta
%
% the vector r = Ap* + c  makes the constraints  Ap + c = r compatible 
% with the trust region for the original problem
%
%              minimize     1/2 p'* G * p  + g'*p  
%             subject to    Ap + c = r,   ||p|| <= Delta
%
% The computation is performed by means of Powell's dogleg method. Hence,
% two calculatons are needed: the Cauchy point and a Newton-like direction
% given by 
% 
%                 p^U = argmin || Ap + c ||^2
%
% i.e. the unconstrained minimizer of || Ap + c ||^2. Since the Hessian 
% of the quadratic is convex; we choose the minimum norm solution.
%
% Powell's method explores the path formed by the Cauchy point and the
% Newton-like direction. 
%
%    Jose Luis Morales
%    Department of Mathematics
%    ITAM,  2015
%
%--------------------------------------------------------------------------

Delta = Delta*0.8;
%
[ m, n ] = size(A);
%
% ... Compute Newton-like step 
%
z = zeros(m,1);
K = [ speye(n)       A'             ; 
        A   spdiags( z, 0:0, m, m) ];
%
[ L, D, P, S ]  = ldl(K);
%
p_B = MinNorm ( c, m, n, L, D, P, S );

if norm(p_B) <= Delta                  % exit: Newton step is inside the TR
    p = p_B;
    type = 'B';
else                                   % form the dogleg path
    B = A'*A;
    g = A'*c;
    p_U = - g*(g'*g)/(g'*B*g);   
    
    %  set the Cauchy point  
    if norm(p_U) > Delta             
        p = (p_U/norm(p_U))*Delta;     
        type = 'C';
     else
         d = p_B - p_U;
         a = d'*d;                     % take a more profitable dogleg step
         b = 2*p_U'*d;
         c = p_U'*p_U - Delta^2;
         
         x = (- b + sqrt( b^2 - 4*a*c))/(2*a);
         p =  p_U + x*d;
         type = 'D';
     end
end
%  
%--------------------------------------------------------------------------
%
function z = MinNorm ( r, m, n, L, D, P, S );
%
% ... This routine computes the minimum norm solution of the system
%
%                      Ax + c = 0
% 
% by solving the system 
%
%                   | I   A' | | x |      | 0 |
%                   | ------ | | - |   =  |---|
%                   | A   0  | | d |      |-c |
%
%--------------------------------------------------------------------------

rhs = [ zeros(n,1) ; r ];

z = P'*(S*rhs);
z = L\z;
z = D\z;
z = L'\z;
z = P'\z;
z = S*z;
%           resize z
z = z(1:n);









