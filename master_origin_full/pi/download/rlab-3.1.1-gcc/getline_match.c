/* getline_match.c */

/* Source code for the RLaB string manipulating function.
 *
 * SEE THE RLABPLUS MANUAL FOR UP-TO-DATE DESCRITPION
 */

/*  This file is a part of RLaB+rlabplus ("Our"-LaB)
    Copyright (C) 2015 M. Kostrun

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    See the file ./COPYING
 ***********************************************************************/

// **************************************************
//
//
// processing of match-compatible patterns
//
//
// **************************************************
char * process_rlab_string_pattern(char * str)
{
  char *pattern=string_buff;
  char *processed_pattern=0;
  int   len=isvalidstring(str);

  if (len<1)
    return 0;

  if (str[0]=='\'')
  {
    if (len > 1)
    {
      // if first character in delimiter is
      //  ' -> this indicates match-compliant pattern, pass the rest of the string
      // as match-compliant expression;
      // otherwise, assume it is a plain string, so pass it embraced in ''
      pattern = (char *) GC_MALLOC (len * sizeof (char));
      memcpy (pattern, str+1, len);
    }
    else
    {
      // user has provided single ' as pattern: convert it so
      // that Pattern matching can work with it as a simple string
      pattern = (char *) GC_MALLOC (3 * sizeof (char));
      pattern[0]='\"';
      pattern[1]='\'';
      pattern[2]='\"';
      pattern[3]='\0';
    }
  }
  else
  {
    pattern = (char *) GC_MALLOC ((len + 3) * sizeof (char));
    pattern[0]='\'';
    memcpy (pattern+1, str, len);
    pattern[len+1]='\'';
    pattern[len+2]='\0';
  }

  processed_pattern = (char *) patmaker(pattern);

  if (MatchError != MATCH_SUCCESS)
  {
    warning_1("Unable to process the pattern."
        " Please read the manual while I stop for a moment!\n");

    if (processed_pattern)
      GC_FREE(processed_pattern);

    processed_pattern = 0;
  }

  return processed_pattern;
}


// **************************************************
//
//
// strsplt related functions
//
//
// **************************************************

static MDS * string_split_string (char * str, int len, char * dlm)
{
  MDS * rval=0;
  int   i, dlmlen=isvalidstring(dlm);
  char *pattern=0;
  const char *processed_pattern=0;

  if (dlmlen<1)
    return mds_CreateScalar(str);

  // if first character in delimiter is
  //  ' -> this indicates match-compliant pattern, pass the rest of the string as a pattern
  // otherwise assume it is a plain string, so pass it embraced in ''
  if ((dlm[0]=='\'') && (dlmlen>1))
  {
    pattern = (char *) GC_MALLOC (dlmlen * sizeof (char));
    for (i=0;i<dlmlen-1;i++)
      pattern[i] = dlm[i+1];
    pattern[i]='\0';
  }
  else
  {
//     fprintf(stderr, "string_split_string: 0\n");
    pattern = (char *) GC_MALLOC ((dlmlen + 3) * sizeof (char));
    pattern[0]='\"';
    for (i=0;i<dlmlen;i++)
      pattern[i+1] = dlm[i];
    pattern[dlmlen+1]='\"';
    pattern[dlmlen+2]='\0';
  }

//   fprintf(stderr, "string_split_string: 1\n");

  processed_pattern = patmaker(pattern);

//   fprintf(stderr, "string_split_string: 2\n");
  if ( MatchError != MATCH_SUCCESS)
  {
    warning_1("Unable to process the pattern."
        " Please read the manual while I stop for a moment!\n");
    goto _exit_string_split_mds;
  }

  i=0;
  int ilast=0,icol=0,ncols = 20;
  rval=mds_Create(1,ncols);
  while(i<len)
  {
//     fprintf(stderr, "string_split_string: 3: %i\n", i);

    int iptr = 0;
    char *w = (char *) str+i;

//     fprintf(stderr, "w =%s\n", w);

    if (match (strlen(w), w,  &iptr, processed_pattern) == MATCH_SUCCESS)
    {
//       fprintf(stderr, "string_split_string: 4: i=%i iptr=%i\n", i, iptr);
        // (i):(i+iptr-1) contains delimiter
        // ilast:(i-1) contains the string since last time delimiter was found
      if (icol >= ncols-1)
      {
        ncols += 10;
        mds_Extend (rval, 1, ncols);
      }
      MdsV0(rval,icol) = cpnstr ((char *) str+ilast, (i-ilast));
      icol++;
      i += iptr;  // skip delimiter
      ilast = i;  // remember for the next
      iptr=0;     // look at 'w'
    }
    else
    {
//       fprintf(stderr, "string_split_string: 5: %i\n", i);
      if (i==len-1)
      {
          // we got so far without the match
        MdsV0(rval,icol) = cpnstr ((char *) str+ilast, (len-ilast));
        icol++;
      }
      i++;  // nothing was found starting from 'i' proceed to next.
    }
  }
  mds_Extend (rval, 1, icol);

_exit_string_split_mds:

  // clean-up stack
  freedm ();

  // construct pattern for match
  if (pattern)
    GC_FREE(pattern);

  return rval;
}

static MDS * string_split_mdr (char * str, int len, MDR * fwidths)
{
  MDR *rval=0;
  int nf0=SIZE(fwidths),nf;

  if (nf0<1)
    return mds_CreateScalar(str);

  /* Split based upon field width. */
  int i, j, fw, tlen = 0;

  if (nf0==1)
  {
    /* Split the string into constant length fields. */
    fw = (int) MdrV0(fwidths,0);          /* Field Width */
    nf = (int) ((len-1) / fw) + 1;        /* Number Fields */

    if (fw>len)
      return mds_CreateScalar(str);

    /* Create string matrix to hold split results */
    rval = mds_Create (1, nf);
    j = 0;
    for (i=0;i<nf; i++)
    {
      MdsV0 (rval, i) = cpnstr ((str + j), fw);
      j += fw;
    }
    return rval;
  }

  /* Split the string into variable length fields. */
  tlen=i=0;
  while(i<nf0)
  {
    tlen += mdiV0(fwidths,i);
    if (tlen >= len)
      break;
    i++;
  }
  nf = i+1; // number of fields which total length is equal to or greater then the length of the string

  rval = mds_Create (1, nf);
  j = 0;
  for (i=0; i<nf; i++)
  {
    if (i<nf0)
      fw = mdiV0 (fwidths, i);
    else
      fw = len - j;

    MdsV0 (rval, i) = cpnstr ((str + j), fw);
    j = j + fw;
  }
  return rval;
}


// **************************************************
//
//
// strindex related functions
//
//
// **************************************************

int string_strindex_string(char *str, int len, char *dlm, int * flen)
{
  int rval=0;
  int   i, dlmlen=isvalidstring(dlm);
  char *pattern=0;
  char *processed_pattern=0;

  if (dlmlen<1)
    goto _exit_string_strindex_string;

  // if first character in delimiter is
  //  ' -> this indicates match-compliant pattern, pass the rest of the string as a pattern
  // otherwise assume it is a plain string, so pass it embraced in ''
  if (dlm[0]=='\'')
  {
    pattern = (char *) GC_MALLOC (dlmlen * sizeof (char));
    for (i=0;i<dlmlen-1;i++)
      pattern[i] = dlm[i+1];
    pattern[i]='\0';
  }
  else
  {
    pattern = (char *) GC_MALLOC ((dlmlen + 3) * sizeof (char));
    pattern[0]='\'';
    for (i=0;i<dlmlen;i++)
      pattern[i+1] = dlm[i];
    pattern[dlmlen+1]='\'';
    pattern[dlmlen+2]='\0';
  }

  processed_pattern = (char *) patmaker(pattern);
  if (MatchError != MATCH_SUCCESS)
  {
    warning_1(THIS_SOLVER "Unable to process the pattern."
        " Please read the manual while I stop for a moment!\n");
    goto _exit_string_strindex_string;
  }

  i=0;
  while(i<len)
  {
    int iptr = 0;
    char *w = (char *) str+i;

    i++;
    if (match (strlen(w), w,  &iptr, (const char *)processed_pattern) == MATCH_SUCCESS)
    {
      rval = i;

      if (flen)
        *flen = iptr; // make note of the length of the pattern found

      break;
    }
  }

_exit_string_strindex_string:

  freedm ();

  if (pattern)
    GC_FREE(pattern);
  if (processed_pattern)
    GC_FREE(processed_pattern);

  return rval;
}

// **************************************************
//
//
// strtod related functions
//
//
// **************************************************

//
// given array of processed patterns look for one of their first appearance
// in string
//    if found
//        return its start, and set iptr_ex to its length
//    otherwise
//        return -1
int mds_find_first_left_pattern_in_string(char * str, MDS * processed_pattern, int * iptr_ex)
{
  int k, i, iptr;
  int len=isvalidstring(str);
  int numpat=SIZE(processed_pattern);

  for (i=0; i<len; i++)
  {
    iptr=i;

    for (k=0; k<numpat; k++)
    {
      if (isvalidstring(MdsV0(processed_pattern,k))<1)
        continue;

      if (match (len, str,  &iptr, MdsV0(processed_pattern,k)) == MATCH_SUCCESS)
      {
        if (iptr_ex)
          *iptr_ex = iptr;

        return i;
      }
    }

  } // for (i=0; i<len; i++)

  return -1;
}

int string_find_first_left_pattern_in_string(char * str, char * processed_pattern, int * iptr_ex)
{
  int i=0, len=isvalidstring(str);

  if (len<0)
    return -1;

  if (isvalidstring(processed_pattern)<1)
    return -1;

  while(i<len)
  {
    char *w = (char *) str+i;
    int iptr=0;

    if (match (strlen(w), w,  &iptr, processed_pattern)==MATCH_SUCCESS)
    {
      if (iptr_ex)
        *iptr_ex = iptr;
      return i;
    }

    i++;
  }

  freedm();

  return -1;
}

int processed_pattern_starts_string(char * str, MDS * s)
{
  int i, keep_doing_it, iptr=0;

  int   len=isvalidstring(str);

  if (len<1)
    return iptr;

  do
  {
    keep_doing_it = 0;

    for (i=0; i<SIZE(s); i++)
    {
      char * processed_pattern = MdsV0(s,i);
      if (isvalidstring(processed_pattern)<1)
        continue;
      while (match(len, str,  &iptr, processed_pattern)==MATCH_SUCCESS)
        keep_doing_it = 1;
    }

  }
  while (keep_doing_it);

  freedm();

  return iptr;
}

int does_processed_pattern_end_string(char * str, MDS * s)
{
  int i, keep_doing_it;

  int   wlen=isvalidstring(str);
  int   ileft,iright,iptr=0;

  if (wlen<0)
    return iptr;

  ileft = wlen-1;
  iright = wlen;

  do
  {
    keep_doing_it = 0;

    for (i=0; i<SIZE(s); i++)
    {

      char * processed_pattern = MdsV0(s,i);
      if (isvalidstring(processed_pattern)<1)
        continue;

      while (ileft>=0)
      {
        iptr = ileft;
        if (match(iright, str,  &iptr, processed_pattern)==MATCH_SUCCESS)
        {
          if (iptr == iright)
          {
            iright = ileft;
            keep_doing_it=1;
          }
        }
        ileft--;
      }
    }
  }
  while (keep_doing_it && ileft>=0);

  freedm();

  return iright;
}


char * remove_pattern_from_string(MDS * processed_pattern, char * str)
{
  int i_start, i_new_start, i_end, i_clean;
  char *w=0;

  int len=isvalidstring(str);

  if (len<1)
    return 0;

  char *str_clean = (char *) GC_malloc(len+1);

  i_start=0;
  i_clean=0;
  while (i_start<=len)
  {
    w = str + i_start;
    i_end = mds_find_first_left_pattern_in_string(w, processed_pattern, &i_new_start);
    if (i_end == -1)
      break;

    memcpy(str_clean+i_clean, w, i_end);
    i_clean += i_end;
    i_start += i_new_start;
  }
  memcpy(str_clean+i_clean, w, isvalidstring(w));
  *(str_clean + i_clean + isvalidstring(w)) = '\0';

  return str_clean;
}


void strtod_split_string ( MDR **rval, int * rowidx, int * icol, char * str, int len, char * dlm,
                           MDR *set_usecols, int * max_icol )
{
  int   i, dlmlen=isvalidstring(dlm), idx_set_usecols=0;
  char  *end;
  char *processed_pattern=0;

  if (dlmlen<1)
  {
    return;
  }

  processed_pattern = process_rlab_string_pattern(dlm);

  if (!processed_pattern)
    return;

  int ncols = 20;

  *icol = 0;

  if (set_usecols)
  {
    ncols = SIZE(set_usecols);

    // skip column indices below lower end
    while (mdiV0(set_usecols, idx_set_usecols)<0)
      idx_set_usecols++;
  }

  if (!(*rval))
  {
    *rval=mdr_Create(1,ncols);
    *rowidx=0;
  }

  i=0;
  int ilast=0;
  double dummy;
  while(i<len)
  {
    int iptr = 0;
    char *w = (char *) str+i;

    if (match (strlen(w), w,  &iptr, processed_pattern) == MATCH_SUCCESS)
    {
      dummy = create_nan();

      // (i):(i+iptr-1) contains delimiter
      // ilast:(i-1) contains the string since last time delimiter was found
      ncols = MNC(*rval);
      if (set_usecols)
      {
        if (idx_set_usecols >= ncols-1)
        {
          ncols += 10;
          mdr_Extend (*rval, MNR(*rval), ncols);
        }
      }
      else
      {
        if (*icol >= ncols-1)
        {
          ncols += 10;
          mdr_Extend (*rval, MNR(*rval), ncols);
        }
      }

      if (ilast<i)
      {
        // set the beginning of the delimiter to \0
        char c = *w;
        *w = '\0';

        // read the data nonetheless
        if (ilast<i)
        {
          dummy = strtod ((char *) str+ilast, &end);
          if ( (end == str+ilast) && (dummy==0) )
          {
            dummy = create_nan();
          }
        }

        // put w to its original form
        *w = c;
      }

      if (set_usecols)
      {
        // move to next field, only if current column is in the list of columns we want
        if((mdiV0(set_usecols, idx_set_usecols)-1)==*icol)
        {
          Mdr0(*rval,*rowidx, idx_set_usecols) = dummy;
          idx_set_usecols++;

          if (idx_set_usecols==SIZE(set_usecols))
          {
            *icol = idx_set_usecols;
            break;
          }
        }
      }
      else
      {
        Mdr0(*rval,*rowidx, *icol) = dummy;
      }

      (*icol)++;
      i += iptr;  // skip delimiter
      ilast = i;  // remember for the next
    }
    else
    {
      //
      // are we looking at the last string segment in the line?
      //
      if (i==len-1)
      {
        // process it
        dummy = strtod ((char *) str+ilast, &end);
        if (!set_usecols)
        {
          // we got so far without the match
          if ( (end == str+ilast) && (dummy==0) )
          {
            Mdr0(*rval,*rowidx, *icol) = create_nan();
          }
          else
          {
            Mdr0(*rval,*rowidx, *icol) = dummy;
          }
          (*icol)++;
        }
        else
        {
          if((mdiV0(set_usecols, idx_set_usecols)-1)==*icol)
          {
            Mdr0(*rval,*rowidx, idx_set_usecols) = dummy;
            *icol = idx_set_usecols+1;
          }
        }

      } // if ((i==len-1)&&(*icol>0))

      i++;  // nothing was found starting from 'i' proceed to next.
    }
  }

  if ((ilast==i)&&(i==len))
  {
    // data ended with delimiter
    if (set_usecols)
    {
      // move to next field, only if current column is in the list of columns we want
      if((mdiV0(set_usecols, idx_set_usecols)-1)==*icol)
      {
        Mdr0(*rval,*rowidx, idx_set_usecols) = create_nan();
        idx_set_usecols++;

        if (idx_set_usecols==SIZE(set_usecols))
        {
          *icol = idx_set_usecols;
        }
      }
    }
    else
    {
      Mdr0(*rval,*rowidx, *icol) = create_nan();
    }

    (*icol)++;
  }

  if (*icol > *max_icol)
    *max_icol = *icol;

//   fprintf(stderr, "strtod_split_string: *icol = %i\n", *icol);

  if ( (*rowidx==0) || ((*rowidx>0)&&(*max_icol!=MNC(*rval))) )
    mdr_Extend (*rval, MNR(*rval), *max_icol);

  if (processed_pattern)
    GC_FREE(processed_pattern);

  // clean-up stack
  freedm ();

  return;
}

//
// NOTE: string_buff is used to assemble string
//
char * replace_processed_pattern_in_string(char * rep, char * processed_pattern, char * str, int *k)
{
  int i_left=0, len=isvalidstring(str), iptr, i_clean=0, i_off, lrep = isvalidstring(rep);
  static char *str_clean=0;
  char *w=0;

  (*k)=0;

  if (len<1)
    return str_clean;

  str_clean=string_buff;

  while (i_left<=len)
  {
    w = str + i_left;
    i_off = string_find_first_left_pattern_in_string(w, processed_pattern, &iptr);

    if (i_off == -1) // pattern not found. get out of here
      break;

    (*k)++;

    // copy segment of w prior to the pattern
    memcpy(str_clean+i_clean, w, i_off);

    // now instead of the pattern put 'rep' there
    if (lrep>0)
    {
      memcpy(str_clean+i_clean+i_off, rep, lrep);
      i_clean += lrep;
    }

    // move 'left' pointer in 'w' after the found pattern
    i_left  += i_off + iptr;
    i_clean += i_off;
  }
  memcpy(str_clean+i_clean, w, isvalidstring(w));

  *(str_clean + i_clean + isvalidstring(w)) = '\0';

  return cpstr(str_clean);
}

extern void chomp_string(char *s);

#define RLAB_READM_DEFAULT_NRSTART    100
#undef THIS_SOLVER
#define THIS_SOLVER "readm"
MDR * mdr_ReadGeneric (FILE * fn, int block_size,
                       int iskiprows, char *delim, int min_line_len, int join_rows, char * join_csp,
                       MDS *comment, MDS *note, MDS *lstrip, MDS *grep,
                       MDR *userows, MDR *usecols, MDS *start, MDS *stop, int skipdatarows)
{
  int i, j, join_rows_counter=0;

  MDR *rval=0, *set_userows=0, *set_usecols=0;
  MDS *comment_pattern=0, *note_pattern=0, *lstrip_pattern=0, *start_pattern=0, *stop_pattern=0, *grep_pattern=0;

  int start_pattern_notfound=1, stop_pattern_notfound=1, max_icol=0;

  // process strings to create patterns
  if (comment)
  {
    comment_pattern = mds_Create(1,SIZE(comment));
    j=0;
    for (i=0; i<SIZE(comment);i++)
    {
      if (isvalidstring(MdsV0(comment,i))>0)
        MdsV0(comment_pattern,j++) = process_rlab_string_pattern(MdsV0(comment,i));
    }
    if (j)
      mds_Extend(comment_pattern,1,j);
    else
    {
      mds_Destroy(comment_pattern);
      comment_pattern = 0;
    }
  }

  // process strings to create patterns
  if (note)
  {
    note_pattern = mds_Create(1,SIZE(note));
    j=0;
    for (i=0; i<SIZE(note);i++)
    {
      if (isvalidstring(MdsV0(note,i))>0)
        MdsV0(note_pattern,j++) = process_rlab_string_pattern(MdsV0(note,i));
    }
    if (j)
      mds_Extend(note_pattern,1,j);
    else
    {
      mds_Destroy(note_pattern);
      note_pattern = 0;
    }
  }

  // process strings to create patterns
  if (lstrip)
  {
    lstrip_pattern = mds_Create(1,SIZE(lstrip));
    j=0;
    for (i=0; i<SIZE(lstrip);i++)
    {
      if (isvalidstring(MdsV0(lstrip,i))>0)
        MdsV0(lstrip_pattern,j++) = process_rlab_string_pattern(MdsV0(lstrip,i));
    }
    if (j)
      mds_Extend(lstrip_pattern,1,j);
    else
    {
      mds_Destroy(lstrip_pattern);
      lstrip_pattern = 0;
    }
  }

  // process strings to create patterns
  if (start)
  {
    start_pattern = mds_Create(1,SIZE(start));
    j=0;
    for (i=0; i<SIZE(start);i++)
    {
      if (isvalidstring(MdsV0(start,i))>0)
        MdsV0(start_pattern,j++) = process_rlab_string_pattern(MdsV0(start,i));
    }
    if (j)
      mds_Extend(start_pattern,1,j);
    else
    {
      mds_Destroy(start_pattern);
      start_pattern = 0;
    }
  }

  // process strings to create patterns
  if (stop)
  {
    stop_pattern = mds_Create(1,SIZE(stop));
    j=0;
    for (i=0; i<SIZE(stop);i++)
    {
      if (isvalidstring(MdsV0(stop,i))>0)
        MdsV0(stop_pattern,j++) = process_rlab_string_pattern(MdsV0(stop,i));
    }
    if (j)
      mds_Extend(stop_pattern,1,j);
    else
    {
      mds_Destroy(stop_pattern);
      stop_pattern = 0;
    }
  }

  // process strings to create patterns
  if (grep)
  {
    grep_pattern = mds_Create(1,SIZE(grep));
    j=0;
    for (i=0; i<SIZE(grep);i++)
    {
      if (isvalidstring(MdsV0(grep,i))>0)
        MdsV0(grep_pattern,j++) = process_rlab_string_pattern(MdsV0(grep,i));
    }
    if (j)
      mds_Extend(grep_pattern,1,j);
    else
    {
      mds_Destroy(grep_pattern);
      grep_pattern = 0;
    }
  }

  // selected rows:
  if (SIZE(usecols)>0)
    set_usecols = mdr_VectorSet(usecols);

  if (SIZE(set_usecols)<1)
  {
    mdr_Destroy(set_usecols);
    set_usecols=0;
    usecols=0;
  }

  // selected columns
  if (SIZE(userows)>0)
    set_userows = mdr_VectorSet(userows);
  if (SIZE(set_userows)<1)
  {
    mdr_Destroy(set_userows);
    set_userows=0;
    userows=0;
  }

  int nr_start = RLAB_READM_DEFAULT_NRSTART;
  int idx_set_userows=0;

  if (set_userows)
    nr_start = MIN( RLAB_READM_DEFAULT_NRSTART, mdiV1(set_userows,SIZE(set_userows)) );

  if (set_usecols)
  {
    rval = mdr_Create(nr_start,SIZE(set_usecols));
  }
  else
    rval = mdr_Create(nr_start,20);

  mdr_Nan(rval);
  int irow=0,icol=0, ifoundcomment=-1, datarow=0;
  char c=0, *joined_str=0;

  j=-1;
  while (!feof (fn))
  {
    // we read in system wide buffer
    char * str=string_buff;

    // line counter
    j++;

    // read line from file
    if (!fgets (str, MAX_STRING_BUFF, fn))
      break;

    // strip '\n' and '\r' from end of the string
    chomp_string(str);

    // do we skip it?
    if (j<iskiprows)
      continue;

    // do we grep it?
    if (grep_pattern)
    {
      if (mds_find_first_left_pattern_in_string(str, grep_pattern, NULL) < 0)
        continue;
    }

    char * str_clean=0;
    int len_clean;
    int len=isvalidstring(str);

    // did user provide list of rows they wants processed?
    if (set_userows)
    {

      if (j==0)
      {
        while ((j>(mdiV0(set_userows,idx_set_userows)-1)) && (idx_set_userows < SIZE(set_userows)) )
          idx_set_userows++;
      }

      if (idx_set_userows == SIZE(set_userows))
        break;

      if (j!=(mdiV0(set_userows,idx_set_userows)-1))
        continue;

      idx_set_userows++;
    }

    // do we check for minimal length of the string
    if (len < min_line_len)
      continue;

    // if user provided start pattern: then we look for it before we can
    // process the lines
    if (start_pattern && start_pattern_notfound==1)
    {
      if (processed_pattern_starts_string(str, start_pattern) < 1)
        continue;

      start_pattern_notfound=0;
      if (stop_pattern)
        stop_pattern_notfound=1;
      continue;
    }

    // if user provided stop pattern: then we stop processing the lines
    if (stop_pattern && stop_pattern_notfound==1)
    {
      if (processed_pattern_starts_string(str, stop_pattern) > 0)
      {
        if (!start_pattern)
          break;

        stop_pattern_notfound=0;
        start_pattern_notfound=1;
      }
    }

    if (!stop_pattern_notfound)
      continue;

    // did we implement comments?
    //    if so, temporarily replace the first character of comment with '\n'
    //    to prevent match processing past it
    if (comment_pattern)
    {
      ifoundcomment = mds_find_first_left_pattern_in_string(str,comment_pattern,NULL);
      if (ifoundcomment >= 0)
      {
          // we have found comment
        if (ifoundcomment < min_line_len)
        {
            //  what is left after the comment is removed is not long enough
            //  to warrant further processing
          continue;
        }

          // we have modified string; update its length
        c = str[ifoundcomment];
        str[ifoundcomment] = '\0';
        len = isvalidstring(str);
      }
    }

    //
    // at this point we have a row for sure
    //
    datarow++; // count it

    // user may want to skip first few data points for what ever reason:
    //  we do it here to avoid rewriting large arrays of data on small processors
    if (skipdatarows)
      if (datarow <= skipdatarows)
        continue;

    if (join_rows > 1)
    {
      join_rows_counter++;
      if (join_rows_counter > 1)
        string_concat(&joined_str, join_csp);
      string_concat(&joined_str, str);
      if (join_rows_counter < join_rows)
        continue;

        // fix 'str' before replacing it with another pointer
      if (ifoundcomment > -1)
      {
        str[ifoundcomment] = c;
        ifoundcomment = -1;
      }
      str = joined_str;
      join_rows_counter=0;
      len = isvalidstring(str);
    }

    //
    // processing depends whether the notes are provided
    //
    if (note_pattern)
    {
      // notes are removed from string before the strings are processed
      str_clean = remove_pattern_from_string(note_pattern, str);
      len_clean = isvalidstring(str_clean);

      if (len_clean < min_line_len)
      {
        GC_FREE (str_clean);
        continue;
      }

      // fix 'str' before replacing it with another pointer
      if (ifoundcomment > -1)
      {
        str[ifoundcomment] = c;
        ifoundcomment = -1;
      }

      str = str_clean;
      len = len_clean;
    }

    // finally: apply lstrip to remove anything superficial in the
    // string prior to processing
    int iblank = 0;
    if (lstrip_pattern)
    {
      iblank = processed_pattern_starts_string(str, lstrip_pattern);
    }

    // process string
    strtod_split_string (&rval, &irow, &icol, str+iblank, len-iblank, delim, set_usecols, &max_icol);

    // we implemented comments
    //   revert string to its previous value but only if notes and join_lines have not affected
    //   its position
    if (ifoundcomment > -1)
    {
      str[ifoundcomment] = c;
      ifoundcomment = -1;
    }

    if (str_clean)
    {
      GC_FREE (str_clean);
      str_clean = 0;
    }
    if (joined_str)
    {
      GC_FREE(joined_str);
      joined_str = 0;
    }

    // we check for minimal length of the string
    if (min_line_len)
    {
      // we ignore empty lines
      if (!icol)
        continue;
    }
    irow++; // next row

    // extend matrix with more rows
    if (irow >= MNR(rval) && icol>0)
      mdr_Extend (rval, MNR(rval)+20, icol);
    if (!rval)
      break;
  }

  //
  // adjust size of the return matrix
  //
  if (irow > 0)
    mdr_Extend(rval, irow, MNC(rval));
  else
  {
    mdr_Destroy(rval);
    rval = 0;
  }

  // cleanup
  if (comment_pattern)
    mds_Destroy (comment_pattern);
  if (note_pattern)
    mds_Destroy (note_pattern);
  if (lstrip_pattern)
    mds_Destroy (lstrip_pattern);
  if (set_usecols)
    mdr_Destroy(set_usecols);
  if (set_userows)
    mdr_Destroy(set_userows);
  if (start_pattern)
    mds_Destroy(start_pattern);
  if (stop_pattern)
    mds_Destroy(stop_pattern);
  if (grep_pattern)
    mds_Destroy(grep_pattern);

  return rval;
}















