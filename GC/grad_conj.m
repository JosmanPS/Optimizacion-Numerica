function x = grad_conj( H, g)

% =========================================================
%
% Este programa encuentra la solucion (x) del problema:
%				Hx = g
% a traves del metodo de gradiente conjugado. Notar que este
% problema es equivalente al de encontrar el minimos sobre
% una funcion cuadratica
%
% 8 Marzo 2015
%
% Modificación: 24 de Septiembre del 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% Entrada:
% 	- H : matriz del problema (Hessiana)
% 	- g : vector de respuestas (gradiente)
%
% Salida:
% 	- x : el vector solucion al problema 
%
% =========================================================

% Calculamos los valores iniciales
n = length(g);
maxit = 2*n;
x = zeros(n,1);
r = g - H * x;
d = r;
iter = 0;

% Guardamos los resultados para reducir operaciones
rTr = r' * r;
Hd = H * d;
dHd = d' * Hd;
TOL_r0 = 1e-16 * rTr;
TOL2 = 1e-6;

while(iter < maxit && rTr > TOL_r0 && dHd > TOL2)

	% Actualizamos los valores
	alpha = rTr / dHd;
	x = x + alpha * d;
	r = r - alpha * Hd;
	Beta = (r' * r) / rTr;
	d = r + Beta * d;

	% Guardamos los resultados para reducir operaciones
	rTr = r' * r;
	Hd = H * d;
	dHd = d' * Hd;

	% Siguiente iteracion
	iter = iter + 1;

end