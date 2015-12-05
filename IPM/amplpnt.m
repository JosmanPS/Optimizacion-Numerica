function [x,s,l] = amplpnt;

global gAmplStub x0 y0 s0

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

x = x0;
l = y0;
s = s0;