%  amplstub: evaluates the function, gradients & Hessian
%
% input parameters:	    x = current primal variables
%			            l = current Lagrange multipliers
%
% output parameters:	f = objective value
% (colm vectors)    	g = objective gradient
%			            c = constraint functions
%			            A = Jacobian of constraints
%			            W = Hessian of Lagrangian
%


function [f, cc, g, AA, W] = amplstub(x, s, l);

global gAmplStub irows ibnd signs Dsigns d E compressDual

if( isempty( gAmplStub ) )
    while( 1 ) 
	stub = input( 'Choose an ampl stub file : ', 's' );
	if( exist( stub ) ) 
	    break;
	else
	    fprintf( 'The stub file %s does not seem to exist\n', stub );
	end
    end
    use_ampl_stub( stub );
end

[f, c] = spamfunc(x, 0);

neq = length(l) - length(s);
nineq = length(s);
cc = [zeros(neq,1);s] + signs .* [ c(irows); x(ibnd) ]  -  d;

[g, Jac] = spamfunc(x, 1);
A =  Dsigns * [Jac(irows,:); E];
AA = [A';[sparse(nineq, neq),spdiags(s,0,length(s),length(s))]];

y = -compressDual * l;

W = spamfunc(y);
