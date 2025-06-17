static int  icount=0;

static H5L_iterate_t op_func_h5ls (hid_t loc_id, const char *name, const H5O_info_t *info,
            void *v)
{
  MDS * sval = v;

  if (sval)
    Mds0(sval,icount,0) = cpstr((char *)name);

  icount++;

  return 0;
}

#undef THIS_SOLVER
#define THIS_SOLVER "h5ls"
Ent *
ent_h5ls (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  char *filename=0, *objname=0;

  Rfile *rf=0;

  int close_after_rw=0, ires, i, rec=0, find_content=0;
  H5O_info_t object1_info,object2_info;
  hid_t work_obj_id;

  MDS *sval=0;

  // Check n_args
  if (nargs != 1 && nargs != 2 && nargs != 3)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED "\n");
    goto _exit_h5ls;
  }

  //
  // get file name from argument list, always the 1st argument
  //
  e1 = bltin_get_ent (args[0]);
  filename = class_char_pointer (e1);
  if (isvalidstring(filename)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_h5ls;
  }

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "r" );
    }
    close_after_rw = 1;
  }
  if (!rf)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5ls;
  }

  if (rf->filetype != RFILE_H5)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5ls;
  }


  //
  // get object name
  //
  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if(ent_type(e2)==MATRIX_DENSE_STRING)
    {
      objname = class_char_pointer (e2);
      if (isvalidstring(objname)<1)
      {
        printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
        goto _exit_h5ls;
      }
    }
  }

  if (!objname)
    objname = "/";

  // if object name ends with "/" and is a group, then return
  // list of the objects in the group
  find_content = (objname[strlen(objname)-1] == '/');

  //
  // do the recursive search?
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    rec = ((int) class_double (e3) > 0);
  }

  // Turn off error handling
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // check if object exists in the file
  ires = H5Oget_info_by_name(rf->h5_file, objname, &object1_info, H5P_DEFAULT);

  if (ires<0)
    goto _exit_h5ls;

  if (object1_info.type == H5O_TYPE_GROUP)
  {
    if (find_content)
    {
      // open the group object
      work_obj_id = H5Oopen(rf->h5_file, objname, H5P_DEFAULT);

      icount = 0;

      // get the count
      if (rec)
        ires = H5Lvisit   (work_obj_id, H5_INDEX_NAME, H5_ITER_INC, (H5L_iterate_t) op_func_h5ls, (void*) NULL);
      else
        ires = H5Literate (work_obj_id, H5_INDEX_NAME, H5_ITER_INC, NULL, (H5L_iterate_t) op_func_h5ls, (void*) NULL);

      if (icount)
      {
        // go over the count and figure the names
        sval = mds_Create(icount, 2);
        icount = 0;

        if (rec)
          ires = H5Lvisit   (work_obj_id, H5_INDEX_NAME, H5_ITER_INC, (H5L_iterate_t) op_func_h5ls, (void*) sval);
        else
          ires = H5Literate (work_obj_id, H5_INDEX_NAME, H5_ITER_INC, NULL, (H5L_iterate_t) op_func_h5ls, (void*) sval);

        // go once more and directly query their type (doing it in op_func gives only "GROUP" ?!)
        for (i=0; i<icount; i++)
        {
          ires = H5Oget_info_by_name(work_obj_id, Mds0(sval,i,0), &object2_info, H5P_DEFAULT);
          if (object2_info.type == H5O_TYPE_DATASET)
            Mds0(sval,i,1) = cpstr("DATASET");
          else if (object2_info.type == H5O_TYPE_GROUP)
            Mds0(sval,i,1) = cpstr("GROUP");
          else if (object2_info.type == H5O_TYPE_NAMED_DATATYPE)
            Mds0(sval,i,1) = cpstr("DATATYPE");
          else
            Mds0(sval,i,1) = cpstr("UNKNOWN");

          // lastly - rewrite the object name
          add2strings2last(objname, &Mds0(sval,i,0));
        }
      }
      H5Oclose(work_obj_id);
    }
    else
    {
      // if it were group just list its name and return it: nothing to do!
      sval = mds_Create(1,2);
      MdsV0(sval,0) = cpstr(objname);
      MdsV0(sval,1) = cpstr("GROUP");
    }
  }
  else if (object1_info.type == H5O_TYPE_DATASET)
  {
    // if it were dataset just list its name and return it: nothing to do!
    sval = mds_Create(1,2);
    MdsV0(sval,0) = cpstr(objname);
    MdsV0(sval,1) = cpstr("DATASET");
  }
  else if (object1_info.type == H5O_TYPE_NAMED_DATATYPE)
  {
    sval = mds_Create(1,2);
    MdsV0(sval,0) = cpstr(objname);
    MdsV0(sval,1) = cpstr("DATATYPE");
  }

_exit_h5ls:

  if (close_after_rw)
    rfile_Destroy(filename);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDS (sval);
}

#undef THIS_SOLVER
#define THIS_SOLVER "h5mv"
Ent *
ent_h5mv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  char *filename=0, *srcname=0, *destname=0;

  Rfile *rf=0;

  int close_after_rw=0, ires;
  H5O_info_t object1_info;
  hid_t src_loc_id;

  double rval=1;

  // Check n_args
  if (nargs != 2 && nargs != 3)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED "\n");
    goto _exit_h5mv;
  }

  //
  // get file name from argument list, always the 1st argument
  //
  e1 = bltin_get_ent (args[0]);
  filename = class_char_pointer (e1);
  if (isvalidstring(filename)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_h5mv;
  }

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "a" );
      close_after_rw = 1;
    }
  }
  if (!rf)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5mv;
  }
  if (rf->filetype != RFILE_H5)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5mv;
  }
  src_loc_id = rf->h5_file;

  rval = 2;

  //
  // get source object name
  //
  e2 = bltin_get_ent (args[1]);
  srcname = class_char_pointer (e2);
  if (isvalidstring(srcname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_h5mv;
  }
  if (srcname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5mv;
  }

  // Turn off error handling
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // check that source exists in the file
  ires = H5Oget_info_by_name(rf->h5_file, srcname, &object1_info, H5P_DEFAULT);
  if (ires<0)
  {
    printf (THIS_SOLVER ": Warning: Source object '%s' does not exist. Cannot continue!\n", srcname);
    goto _exit_h5mv;
  }

  if (nargs == 2)
  {
    H5Ldelete(src_loc_id, srcname, H5P_DEFAULT);

    rval = 0;

    goto _exit_h5mv;
  }

  //
  // third argument is the destination object name
  //
  rval = 3;
  e3 = bltin_get_ent (args[2]);
  destname = class_char_pointer (e3);

  if (isvalidstring(destname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    goto _exit_h5mv;
  }
  if (destname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5mv;
  }

  // source and destination are the same: nothing to do
  if (!strcmp(destname,srcname))
    goto _exit_h5mv;

  // check that destination does not exist in the file
  ires = H5Oget_info_by_name(src_loc_id, destname, &object1_info, H5P_DEFAULT);
  if (ires>=0)
  {
    if (object1_info.type != H5O_TYPE_GROUP)
    {
      printf (THIS_SOLVER ": Warning: Destination object '%s' in file %s exists and is not group. Cannot continue!\n",
              destname, filename);
      goto _exit_h5mv;
    }
  }

  rval = 0;

  // create the destination object
  hid_t gcpl = H5Pcreate (H5P_LINK_CREATE);
  ires = H5Pset_create_intermediate_group (gcpl, 1);

  H5Lmove(src_loc_id, srcname, src_loc_id, destname, gcpl, H5P_DEFAULT );

_exit_h5mv:

  if (close_after_rw)
    close_file_ds(filename);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Double (rval);
}


#undef THIS_SOLVER
#define THIS_SOLVER "h5cp"
Ent * ent_h5cp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  char *filename=0, *filename2=0, *srcname=0, *destname=0;

  Rfile *rf=0, *rf2=0;

  int close_after_rw=0, close_after_rw2=0, ires;
  H5O_info_t object1_info;
  hid_t src_loc_id, dest_loc_id;

  double rval=1;

  // Check n_args
  if (nargs != 3 && nargs != 4)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED "\n");
    goto _exit_h5cp;
  }

  //
  // get file name from argument list, always the 1st argument
  //
  e1 = bltin_get_ent (args[0]);
  filename = class_char_pointer (e1);
  if (isvalidstring(filename)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_h5cp;
  }

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "a" );
      close_after_rw = 1;
    }
  }
  if (!rf)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5cp;
  }
  if (rf->filetype != RFILE_H5)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5cp;
  }
  dest_loc_id = src_loc_id  = rf->h5_file;

  rval = 2;

  //
  // get source object name
  //
  e2 = bltin_get_ent (args[1]);
  srcname = class_char_pointer (e2);
  if (isvalidstring(srcname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_h5cp;
  }
  if (srcname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5cp;
  }

  // Turn off error handling
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // check that source exists in the file
  ires = H5Oget_info_by_name(rf->h5_file, srcname, &object1_info, H5P_DEFAULT);
  if (ires<0)
  {
    printf (THIS_SOLVER ": Warning: Source object '%s' does not exist. Cannot continue!\n", srcname);
    goto _exit_h5cp;
  }

  if (nargs == 3)
  {
    //
    // third argument is the destination object name
    //
    rval = 3;

    e3 = bltin_get_ent (args[2]);
    destname = class_char_pointer (e3);
  }
  else
  {
    //
    // third argument is the destination filename
    //
    rval = 3;

    e3 = bltin_get_ent (args[2]);
    filename2 = class_char_pointer (e3);
    if (isvalidstring(filename2)<1)
    {
      printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
      goto _exit_h5cp;
    }

    //
    // Try and open() the file
    //
    rf2 = rfile_find(filename2);
    if (!rf2)
    {
      // go over supported protocols and try to open the file for
      // close-after writing
      if (    !strncmp(filename2, "h5://", 5)
               ||  !strncmp(filename2, "hdf5://", 7)
               ||  (strstr (filename2, ".h5"))
         )
      {
        // is it HDF5?
        if (H5Fis_hdf5(filename2) > 0)
          rf2 = get_rfile_ds (filename2, RFILE_H5, "a" ); // yes: file2 exists, just append to it!
        else
          rf2 = get_rfile_ds (filename2, RFILE_H5, "w" ); // no:  create file2 as write

        close_after_rw2 = 1;
      }
    }
    if (!rf2)
    {
      printf (THIS_SOLVER  ": Cannot open file %s for write\n", filename2);
      goto _exit_h5cp;
    }
    if (rf2->filetype != RFILE_H5)
    {
      printf (THIS_SOLVER  ": Cannot open file %s for write\n", filename2);
      goto _exit_h5cp;
    }
    dest_loc_id =  rf2->h5_file;

    //
    //
    //
    rval = 4;

    e4 = bltin_get_ent (args[3]);
    destname = class_char_pointer (e4);
  }

  if (isvalidstring(destname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    goto _exit_h5cp;
  }
  if (destname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5cp;
  }

  // source and destination are the same: nothing to do
  if (filename2)
  {
    if (!strcmp(destname,srcname) && !strcmp(filename,filename2))
      goto _exit_h5cp;
  }
  else
  {
    if (!strcmp(destname,srcname))
      goto _exit_h5cp;
  }

  // check that destination does not exist in the file
  ires = H5Oget_info_by_name(dest_loc_id, destname, &object1_info, H5P_DEFAULT);
  if (ires>=0)
  {
    if (object1_info.type != H5O_TYPE_GROUP)
    {
      printf (THIS_SOLVER ": Warning: Destination object '%s' in file %s exists and is not group. Cannot continue!\n",
              destname, filename2);
      goto _exit_h5cp;
    }
  }

  rval = 0;

  hid_t ocpl_id = H5Pcreate (H5P_OBJECT_COPY);
  H5Pset_create_intermediate_group (ocpl_id, 1);

  hid_t lcpl_id = H5Pcreate (H5P_LINK_CREATE);
  H5Pset_create_intermediate_group(lcpl_id, 1);

  ires = H5Ocopy (src_loc_id, srcname, dest_loc_id, destname, ocpl_id, lcpl_id );
  if (ires < 0)
    fprintf(stderr, THIS_SOLVER ": Warning: Object copy operation failed with number %i!\n", ires);

  H5Pclose(ocpl_id);
  H5Pclose(lcpl_id);

_exit_h5cp:

  if (close_after_rw)
    close_file_ds(filename);
  if (close_after_rw2)
    close_file_ds(filename2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Double (rval);
}

#undef THIS_SOLVER
#define THIS_SOLVER "h5ln"
Ent * ent_h5ln (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  char *filename=0, *filename2=0, *srcname=0, *destname=0;

  Rfile *rf=0, *rf2=0;

  int close_after_rw=0, close_after_rw2=0, ires;
  H5O_info_t object1_info;
  hid_t src_loc_id, dest_loc_id;

  double rval=1;

  // Check n_args
  if (nargs != 3 && nargs != 4)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED "\n");
    goto _exit_h5ln;
  }

  //
  // get file name from argument list, always the 1st argument
  //
  e1 = bltin_get_ent (args[0]);
  filename = class_char_pointer (e1);
  if (isvalidstring(filename)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_h5ln;
  }

  //
  // Try and open() the file
  //
  rf = rfile_find(filename);
  if (!rf)
  {
    // go over supported protocols and try to open the file for
    // close-after writing
    if (    !strncmp(filename, "h5://", 5)
        ||  !strncmp(filename, "hdf5://", 7)
        ||  (strstr (filename, ".h5"))
       )
    {
      // it is HDF5
      rf = get_rfile_ds (filename, RFILE_H5, "a" );
      close_after_rw = 1;
    }
  }
  if (!rf)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5ln;
  }
  if (rf->filetype != RFILE_H5)
  {
    printf (THIS_SOLVER  ": Cannot open file %s for read\n", filename);
    goto _exit_h5ln;
  }
  dest_loc_id = src_loc_id  = rf->h5_file;
  filename2 = filename;

  rval = 2;

  //
  // get source object name
  //
  e2 = bltin_get_ent (args[1]);
  srcname = class_char_pointer (e2);
  if (isvalidstring(srcname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_h5ln;
  }
  if (srcname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5ln;
  }

  // Turn off error handling
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);

  // check that source exists in the file
  ires = H5Oget_info_by_name(rf->h5_file, srcname, &object1_info, H5P_DEFAULT);
  if (ires<0)
  {
    printf (THIS_SOLVER ": Warning: Source object '%s' does not exist. Cannot continue!\n", srcname);
    goto _exit_h5ln;
  }

  if (nargs == 3)
  {
    //
    // third argument is the destination object name
    //
    rval = 3;

    e3 = bltin_get_ent (args[2]);
    destname = class_char_pointer (e3);
  }
  else
  {
    //
    // third argument is the destination filename
    //
    rval = 3;

    e3 = bltin_get_ent (args[2]);
    filename2 = class_char_pointer (e3);
    if (isvalidstring(filename2)<1)
    {
      printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
      goto _exit_h5ln;
    }

    //
    // Try and open() the file
    //
    rf2 = rfile_find(filename2);
    if (!rf2)
    {
      // go over supported protocols and try to open the file for
      // close-after writing
      if (    !strncmp(filename2, "h5://", 5)
               ||  !strncmp(filename2, "hdf5://", 7)
               ||  (strstr (filename2, ".h5"))
         )
      {
        // is it HDF5?
        if (H5Fis_hdf5(filename2) > 0)
          rf2 = get_rfile_ds (filename2, RFILE_H5, "a" ); // yes: file2 exists, just append to it!
        else
          rf2 = get_rfile_ds (filename2, RFILE_H5, "w" ); // no:  create file2 as write

        close_after_rw2 = 1;
      }
    }
    if (!rf2)
    {
      printf (THIS_SOLVER  ": Cannot open file %s for write\n", filename2);
      goto _exit_h5ln;
    }
    if (rf2->filetype != RFILE_H5)
    {
      printf (THIS_SOLVER  ": Cannot open file %s for write\n", filename2);
      goto _exit_h5ln;
    }
    dest_loc_id =  rf2->h5_file;

    //
    //
    //
    rval = 4;

    e4 = bltin_get_ent (args[3]);
    destname = class_char_pointer (e4);
  }

  if (isvalidstring(destname)<2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    goto _exit_h5ln;
  }
  if (destname[0]!='/')
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_SCALAR "\n");
    printf (THIS_SOLVER ": Only absolute paths are accepted for HDF5-objects!\n");
    goto _exit_h5ln;
  }

  // source and destination are the same: nothing to do
  if (filename2)
  {
    if (!strcmp(destname,srcname) && !strcmp(filename,filename2))
      goto _exit_h5ln;
  }
  else
  {
    if (!strcmp(destname,srcname))
      goto _exit_h5ln;
  }

  // check that destination does not exist in the file
  ires = H5Oget_info_by_name(dest_loc_id, destname, &object1_info, H5P_DEFAULT);
  if (ires>=0)
  {
    if (object1_info.type != H5O_TYPE_GROUP)
    {
      printf (THIS_SOLVER ": Warning: Destination object '%s' in file %s exists and is not group. Cannot continue!\n",
              destname, filename2);
      goto _exit_h5ln;
    }
  }

  rval = 0;
  hid_t lcpl_id = H5Pcreate (H5P_LINK_CREATE);
  H5Pset_create_intermediate_group (lcpl_id, 1);

  if (nargs==3 || !strcmp(filename,filename2))
  {
    ires = H5Lcreate_soft(srcname, dest_loc_id, destname, lcpl_id, H5P_DEFAULT );
    if (ires < 0)
      fprintf(stderr, THIS_SOLVER ": Warning: creation of soft link failed with number %i!\n", ires);
  }
  else
  {
    ires = H5Lcreate_external(filename, srcname, dest_loc_id, destname, lcpl_id, H5P_DEFAULT );
    if (ires < 0)
      fprintf(stderr, THIS_SOLVER ": Warning: creation of external soft link failed with number %i!\n", ires);
  }

  H5Pclose(lcpl_id);

_exit_h5ln:

  if (close_after_rw)
    close_file_ds(filename);
  if (close_after_rw2)
    close_file_ds(filename2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Double (rval);
}
