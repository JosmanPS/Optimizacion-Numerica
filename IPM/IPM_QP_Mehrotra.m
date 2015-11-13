function [ x, lm, s ] = IPM_QP_Mehrotra(Q, c, A, b, tol, maxiter)

	% ---------------------------------------------------------------
	%
	%   Interior Point Methods (Mehrotra) for Quadratic Programming
	%
	% Solves the problem:
	%        min        c'*x + 1/2*x'*Q*x
	%        s.t.          A*x - b = 0
	%
	% ---------------------------------------------------------------

	%
	% Initial point
	%
	[m, n] = size(A);
	x = ones(n, 1);
	s = x;
	lm = A' \ (Q*x - s + c);
	mu = (x' * s) / n;
	d_gap = mu;
	F1 = Q*x - A'*lm - s + c;
	F2 = A*x - b;
	F3 = x .* s;
	obj = 0.5 * x' * Q * x + c' * x;
	iter = 0;

	%
	% Initial print
	%
	fprintf('iter   ||F1||      ||F2||      ||F3||         OBJ            mu        sigma       alpha \n');
	fprintf('------------------------------------------------------------------------------------------- \n');
	fprintf('%3i  %1.4e  %1.4e  %1.4e  %1.7e  %1.4e \n', iter, norm(F1), norm(F2), norm(F3), obj, mu);

	while d_gap > tol && iter < maxiter

		%
		% Predictor system
		%
		X_inv = sparse(diag(1./x));
		X_inv_S = sparse(diag((1./x) .* s));
		KKT = [ Q + X_inv_S ,         A'     ;
		             A      ,   sparse(m, m)];
		F = [F1 + s; F2];

		%
		% Solve predictor system
		%
		d = backsolve(KKT, -F);
		dx = d(1:n);
		dlm = -d((n+1):(n+m));
		ds = - X_inv_S * dx - s;

		alpha_x = step(x, dx, 1);
		alpha_s = step(s, ds, 1);

		%
		% Compute the central parameter
		%
		mu_aff = (x + alpha_x * dx)' * (s + alpha_s * ds) / n;
		sigma = (mu_aff / mu)^3;

		%
		% Corrector system
		%
		X_inv_F3_c = s - sigma * mu * 1./x;
		F = [F1 + X_inv_F3_c; F2];
		d = backsolve(KKT, -F);
		dx = d(1:n);
		dlm = -d((n+1):(n+m));
		ds = - X_inv_S * dx - X_inv_F3_c;

		alpha_x = step(x, dx, 0.9995);
		alpha_s = step(s, ds, 0.9995);
		alpha = min(alpha_x, alpha_s);

		%
		% Update
		%
		x = x + alpha_x * dx;
		lm = lm + alpha * dlm;
		s = s + alpha_s * ds;

		mu = (x' * s) / n;
		d_gap = mu;
		F1 = Q*x - A'*lm - s + c;
		F2 = A*x - b;
		F3 = x .* s;
		obj = 0.5 * x' * Q * x + c' * x;

		iter = iter + 1;

		%
		% Print
		%
		fprintf('%3i  %1.4e  %1.4e  %1.4e  %1.7e  %1.4e  %1.4e  %1.4e \n', iter, norm(F1), norm(F2), norm(F3), obj, mu, sigma, alpha);

	end

end