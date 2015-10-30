
%
% Esta rutina corre los 59 problemas de nuestro conjunto de prueba para los
% métodos de PCS, tanto el que toma direcciones de Newton completas, como de
% aquel que hace uso de búsqueda lineal sobre la función de mérito
%

% Tomamos la lista de problemas
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
    'byrdsphr';
    'catena';
    'catenary';
    'dixchlng';
    'dtoc1l';
    'dtoc1na';
    'dtoc1nb';
    'dtoc1nc';
    'dtoc1nd';
    'dtoc2';
    'dtoc4';
    'dtoc5';
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
    'hs050';
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
    'orthrdm2';
    'orthrds2';
    'orthrega';
    'orthregb';
    'orthregc';
    'orthregd';
    'orthrgdm';
    'orthrgds'
       };

n_problems = length(list);

% Creamos vectores para guardar los resultados de cada corrida
% Las columnas representan:
%    Tuvo éxito |  nevalf  | CPU time  |  spd  |  iter
newton_results = zeros(n_problems, 5);
newton_results(:, 1) = 1;
merit_results = newton_results;

% Definimos los parámetros
itermax = 50;
tol = 1e-6;
verbose = 0;

%
% Resultados con direcciones de Newton
%
for i = 1:n_problems
    try
        [ ~, ~, ~, ~, iter, feval, time, spd ] = pcs_newton(list{i}, itermax, tol, ...
                                                         verbose);
        newton_results(i, 2) = feval;
        newton_results(i, 3) = time;
        newton_results(i, 4) = spd;
        newton_results(i, 5) = iter;
    catch
        % No se pudo resolver el problema (por recorte)
        newton_results(i, 1) = 0;
    end
end

%
% Resultados con función de mérito
%
for i = 1:n_problems
    try
        [ ~, ~, ~, ~, iter, feval, time, spd ] = pcs_merit(list{i}, itermax, tol, ...
                                                        verbose);
        merit_results(i, 2) = feval;
        merit_results(i, 3) = time;
        merit_results(i, 4) = spd;
        merit_results(i, 5) = iter;
    catch
        % No se pudo resolver el problema
        merit_results(i, 1) = 0;
    end
end

% Buscamos los máximos en feval y time para crear las escalas
max_feval_newton = max(newton_results(:, 2));
max_feval_merit = max(merit_results(:, 2));
max_feval = max(max_feval_newton, max_feval_merit);

max_time_newton = max(newton_results(:, 3));
max_time_merit = max(merit_results(:, 3));
max_time = max(max_time_newton, max_time_merit);

% Ponemos en infinito los valores de los problemas que no resolvimos
index1 = find(newton_results(:, 4) == 0);
index1 = union(find(newton_results(:, 1) == 0), index1);
index1 = union(find(newton_results(:, 5) == itermax), index1);
newton_results(index1, 2:3) = inf;

index2 = find(merit_results(:, 4) == 0);
index2 = union(find(merit_results(:, 1) == 0), index2);
index2 = union(find(merit_results(:, 5) == itermax), index2);
merit_results(index2, 2:3) = inf;

% Definimos el núero de cortes para la gráfica
n_cuts = n_problems;

%
% Gráfica de tiempo de CPU
%
tau = 1:n_cuts;
tau = tau / n_cuts;
tau = tau * max_time;
newton_solved = zeros(1, length(tau));
merit_solved = newton_solved;
for i = 1:n_cuts
    t = tau(i);
    newton_solved(i) = sum(newton_results(:, 3) <= t);
    merit_solved(i) = sum(merit_results(:, 3) <= t);
end

newton_solved = newton_solved / n_problems;
merit_solved = merit_solved / n_problems;

% plot
figure(1);
p = plot(tau, newton_solved, '--', tau, merit_solved, '--');
ylabel('Porcentaje de problemas resueltos');
xlabel('CPU (time)');
p = legend(p, 'PCS-Newton', 'PCS-Mérito');
set(p, 'Location', 'southeast')
        

%
% Gráfica feval
%
tau = 1:n_cuts;
tau = tau / n_cuts;
tau = tau * max_feval;
newton_solved = zeros(1, length(tau));
merit_solved = newton_solved;
for i = 1:n_cuts
    t = tau(i);
    newton_solved(i) = sum(newton_results(:, 2) <= t);
    merit_solved(i) = sum(merit_results(:, 2) <= t);
end

newton_solved = newton_solved / n_problems;
merit_solved = merit_solved / n_problems;

% plot
figure(2);
p = plot(tau, newton_solved, '--', tau, merit_solved, '--');
ylabel('Porcentaje de problemas resueltos');
xlabel('Evaluaciones de la función objetivo');
p = legend(p, 'PCS-Newton', 'PCS-Mérito');
set(p, 'Location', 'southeast')


%
% Gráfico log(time)
%
aux_index = union(index1, index2);
index = 1:n_problems;
index = setdiff(index, aux_index);

log_time = (1 + newton_results(index, 3)) ./ (1 + merit_results(index, 3));
log_time = - log(log_time);
xx = 1:length(log_time);

% plot
figure(3);
bar(xx, log_time);
ylabel('-log( Newton_{CPU} / Mérito_{CPU} )');
