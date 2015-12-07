

list = {
    'airport';
    'aljazzaf';
    'alsotame';
    'batch';
    'bigbank';
    'cantilvr';
    'cb2';
    'cb3';
    'chaconn1';
    'chaconn2';
    'concon';
    'congigmz';
    'csfi2';
    'dipigri';
    'dnieper';
    'eg3';
    'eigmaxb';
    'eigminb';
    'expfitb';
    'gausselm';
    'haifas';
    'hanging';
    'himmelp6';
    'hong';
    'hubfit';
    'loadbal';
    'madsschj';
    'makela2';
    'makela3';
    'mconcon';
    'mifflin1';
    'mifflin2';
    'minmaxbd';
    'minperm';
    'odfits';
    'optcntrl';
    'optprloc';
    % 'polak1';
    'polak3';
    'polak6';
    'prodpl0';
    'prodpl1';
    'rosenmmx';
    'snake';
    'stancmin';
    'swopf';
    'synthes1';
    'trimloss';
    'twobars';
    'zecevic3';
    'zecevic4';
    % 'chenhark';
    % 'clnlbeam';
    % 'cvxbqp1';
    % 'explin';
    % 'explin2';
    % 'expquad';
    % 'jnlbrng1';
    % 'jnlbrng2';
    % 'jnlbrnga';
    % 'jnlbrngb';
    % 'mccormck';
    % 'ncvxbqp1';
    % 'ncvxbqp2';
    % 'ncvxbqp3';
    % 'nobndtor';
    % 'nonscomp';
    % 'obstclae';
    % 'obstclbm';
    % 'pentdi';
    % 'qudlin';
    % 'reading1';
    % 'sineali'
};

TOL = 1e-8;
maxiter = 100;

for i = 1:length(list)
    element = char(list(i));
    [n_x, m, in_iter, f, infeas, spd] = ipm_nl(element, TOL, maxiter);

    if not(spd)
        infeas = 99999;
    end

    fprintf('%s & %5i & %5i & %5i & %1.14e & %1.5e \\\\ \n', ...
            element, n_x, m, in_iter, f, infeas);

end
