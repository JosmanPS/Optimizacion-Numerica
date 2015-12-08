
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

TOL = 1e-8;
maxiter = 200;

for i = 1:length(list)
    element = char(list(i));
    [n, m, iter, f, norm_gL, norm_c] = PCS_RC(element, TOL, maxiter, 0);

    fprintf('%s & %5i & %5i & %5i & %1.14e & %1.5e & %1.5e \\\\ \n', ...
            element, n, m, iter, f, norm_gL, norm_c);
end
