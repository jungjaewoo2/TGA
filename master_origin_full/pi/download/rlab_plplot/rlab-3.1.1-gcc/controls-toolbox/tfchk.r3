//----------------------------------------------------------------------------
//
// tfchk
//
// Syntax: a=tfchk(num,den) or </denc;numc/>=tfchk(num,den)
//
// This routine checks the transfer function given in (num,den) where
// num is the numerator and den is the denominator. The routine passes
// back the transfer function in a list.
//
// If the transfer function is not proper given by (num,den) then
// it returns lenth(numc) = length(denc).
//
// The returned list for a=tfchk(num,den) is:
//
//       a.numc = correct numerator
//       a.denc = correct denominator
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 931016
//----------------------------------------------------------------------------
require isempty

tfchk = function(num,den)
{

   // Check for empty matrices
   if ( isempty(num) ) {
      printf("TFCHK: Warning: Transfer function numerator polynomial is empty.\n");
   }
   if ( isempty(den) ) {
      printf("TFCHK: Warning: Transfer function denominator polynomials is empty.\n");
   }

   // Make sure den is a row vector, den is assumed to be in rows.
   if ( !( (den.nr == 1) || (den.nc == 1)) ) {
      error("TFCHK: Denominator must be a row vector.");
   }

   if ( (den.nr != 1) && (den.nc == 1)) {
      error("TFCHK: Tenominator must be a row vector.");
   }

   if (num.nc > den.nc) {
      error("TFCHK: Transfer function is not proper.");
   }

   // Make num and den lengths equal.
   numc = [zeros(num.nr,den.nc-num.nc),num];
   denc = den;

   return << numc=numc; denc=denc >>;
};

