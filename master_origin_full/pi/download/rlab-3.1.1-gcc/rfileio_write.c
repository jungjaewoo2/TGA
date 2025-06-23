//
// Write out a matrix in a generic sort of format
//
static void mdr_WriteGeneric (MDR * m, Rfile * rf)
{
  int i, j, fwidth, fprec;

  char *eol = rf->eol;
  char *csp = rf->csp;
  char *nan = rf->nan;
  char *inf_pos = rf->inf_pos;
  char *inf_neg = rf->inf_neg;

  MDS  *fmt = rf->fmt;
  FILE *fn  = rf->fileds_f;

  if (SIZE(m)<1)
    return;
  if (!rf)
    return;

  fwidth = get_fwidth ();
  fprec = get_fprec ();

  if(m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < MNR (m); i++)
    {
      for (j=0; j < MNC (m); j++)
      {
        if (j>0)
          fprintf (fn, csp);

        if (isvalidstring(mdsV0_safe(fmt,j))>1)
          fprintf (fn, mdsV0_safe(fmt,j), Mdi0 (m, i, j));
        else
          fprintf (fn, "%*i", fwidth, Mdi0 (m, i, j));
      }
      fprintf (fn, eol);
    }
  }
  else
  {
    double pinf=create_inf(), ninf=-create_inf();
    for (i=0; i<MNR (m); i++)
    {
      for (j=0; j < MNC (m); j++)
      {
        if (j>0)
          fprintf (fn, csp);

        if (isnand(Mdr0 (m, i, j)))
        {
          fprintf(fn, nan);
          continue;
        }

        if (Mdr0 (m, i, j) == pinf)
        {
          fprintf(fn, inf_pos);
          continue;
        }

        if (Mdr0 (m, i, j) == ninf)
        {
          fprintf(fn, inf_neg);
          continue;
        }

        if (isvalidstring(mdsV0_safe(fmt,j))>1)
          fprintf (fn, mdsV0_safe(fmt,j), Mdr0 (m, i, j));
        else
          fprintf (fn, "%*.*f", fwidth, fprec, Mdr0 (m, i, j));
      }
      fprintf (fn, eol);
    }
  }
}

/*
 * Write out a matrix in a generic sort of format
 */

static void mds_WriteGeneric (MDS * m, Rfile * rf)
{
  int i, j, fwidth;

  char *eol = rf->eol;
  char *csp = rf->csp;
  MDS  *fmt = rf->fmt;
  FILE *fn  = rf->fileds_f;

  if (SIZE(m)<1)
    return;
  if (!rf)
    return;

  fwidth = get_fwidth ();

  for (i=0; i < MNR(m); i++)
  {
    for (j=0; j < MNC(m); j++)
    {
      if (j>0)
        fprintf (fn, csp);

      if (isvalidstring(mdsV0_safe(fmt,j))>1)
        fprintf (fn, mdsV0_safe(fmt,j), Mds0 (m, i, j));
      else
        fprintf (fn, "%*s", fwidth, Mds0 (m, i, j));
    }
    fprintf (fn, eol);
  }
}

