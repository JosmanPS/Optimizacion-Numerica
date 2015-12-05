function stub = use_ampl_stub( stub_in )

global gAmplStub irows ibnd Dsigns signs d E compressDual 
global x0 y0 s0

gAmplStub = stub_in;
stub = gAmplStub;

if( ~exist(gAmplStub) )
    warning( sprintf( 'The file %s does not seem to exist.\n', gAmplStub));
end

%[x, xlow, xupp, y, clow, cupp] = spamfunc( gAmplStub ); see corr by j-l
[x, xlow, xupp, y, clow, cupp, ccomp] = spamfunc( gAmplStub );

ieq    = find( clow == cupp );                neq    = length( ieq );
iclow  = find( clow > -inf & clow ~= cupp );  nclow  = length( iclow );
icupp  = find( cupp < inf  & clow ~= cupp );  ncupp  = length( icupp );
% added by j-l to incorporate complementarity problems
icomp  = find(ccomp);                         ncomp  = length( icomp );
ngen   = neq + nclow + ncupp + ncomp;

% assume no fixed variables.
ixlow  = find( xlow > -inf );  nxlow = length( ixlow );
ixupp  = find( xupp < inf );   nxupp = length( ixupp );

II = speye( length( x ) );
E  = [II( ixlow, : ); II( ixupp, :)];
% see correction by j-l
irows = [ieq; iclow; icupp; icomp];
ibnd  = [ixlow; ixupp ];

signs = [ -ones( neq + nclow, 1); ones( ncupp, 1);...
    -ones( nxlow, 1); ones(nxupp, 1) ];
Dsigns = spdiags(signs,0,length(signs), length(signs) );

d = [ clow( ieq ); clow( iclow ); cupp( icupp );...
    xlow( ixlow ); xupp( ixupp ) ] .* signs;

II = speye( length( clow ) );
compressDual = [ II(ieq,:); II(iclow,:); -II(icupp,:);... 
    sparse( nxlow + nxupp, length(clow)) ]';

s0 = ones( nclow + ncupp + nxlow + nxupp, 1 );
x0 = x;
y0 = [y(ieq); ones(size(y(iclow))); ones(size(y(icupp))); ones(nxlow,1); ones(nxupp,1)];
%y0 = max( 1, y0 );


