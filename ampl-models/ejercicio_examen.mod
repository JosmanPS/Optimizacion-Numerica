# AMPL Model by Josman Proudinat
#
# Copyright (C) 2015 ITAM
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.

var x{i in 1..2} := 3-i;

minimize f:
    -2*x[1]+x[2];

subject to cons1:
    (1-x[1])^3+x[2]=0;


display f;
display x;
