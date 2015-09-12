# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: problem 11 in
#   P.T. Boggs and J.W. Tolle,
#   "A strategy for global convergence in a sequential 
#    quadratic programming algorithm",
#   SINUM 26(3), pp. 600-623, 1989.

#   The problem is not convex.

#   SIF input: Ph. Toint, June 1993.

#   classification OOR2-AY-5-3

var x{1..5} := 2.0;

minimize f:
	(x[1]-1.0)^2 + (x[1]-x[2])^2 + (x[2]-x[3])^2 + (x[3]-x[4])^4 + (x[4]-x[5])^4;
subject to cons1:
	x[1]+x[2]^2+x[3]^3 = -2+sqrt(18.0);
subject to cons2:
	x[2]+x[4]-x[3]^2 = -2+sqrt(8.0);
subject to cons3:
	x[1]-x[5] = 2.0;

 
display f;
display x;
