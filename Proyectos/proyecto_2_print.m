
% ----------------------------------------------------------------------
%
% Esta función imprime los resultados del problema evaluado con
% Programación Cuadrática Sucesiva.
%
% Josman, 2015
%
% ----------------------------------------------------------------------

list = {
    'bt11';
    'bt12';
    'bt1';
    'bt2';
    'bt4';
    'bt5';
    'bt6';
    'bt7';
    'bt8';
    'bt9';
    'catena';
    'catenary';
    'dixchlng';
    'dtoc1nb';
    'dtoc1nc';
    'dtoc1nd';
    'dtoc6';
    'eigena2';
    'eigenaco';
    'eigenb2';
    'eigenbco';
    'eigenc2';
    'eigencco';
    'gilbert';
    'hs006';
    'hs007';
    'hs009';
    'hs026';
    'hs027';
    'hs039';
    'hs040';
    'hs046';
    'hs047';
    'hs049';
    'hs061';
    'hs077';
    'hs078';
    'hs079';
    'hs100lnp';
    'hs111lnp';
    'lch';
    'maratos2';
    'maratos';
    'mwright';
    'orthregb';
    'orthrgds'
       };

for i = 1:length(list)
    element = char(list(i));
    [ x, n, m, f, iter, feval, time, spd ] = pcs_newton(element, 50, 1e-6, 0);

    fprintf(['%s  &  %5i  &  %5i  &  %1.14e  &  %3i  &  %3i  &  %1.2e & %1i \\\\ \n' ...
             ''], element, n, m, f, iter, feval, time, spd);
end
