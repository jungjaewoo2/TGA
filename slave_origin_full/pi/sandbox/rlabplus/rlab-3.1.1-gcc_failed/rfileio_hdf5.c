// rfileio_hdf5.c
//
// collection of functions that provide HDF5 functionality within rlab
//
/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 2010  Marijan Kostrun, part of rfileio.c by Ian Searle,

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
   **********************************************************************
*/


void h5_create_groups_for_writing_atomic_object (Rfile * rf, char * objname)
{
  // no error messages:
  if(!RLAB_HDF5_DEBUG)
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  if (rf->h5_flags != H5F_ACC_RDONLY)
  {
    // Requested object ID does not exist, but the file is write/append
    // thus in anticipation of writing create the groups leading to the
    // named object. Use
    //   obj_name = "{/}name0/name1/.../nameN
    // to create groups /name0, /name0/name1 ... to /name0/.../name{N-1} as needed

    // count '/' from the second character
    char *cntsl = objname;
    int  icnt = 1;
    while (cntsl[icnt] != '\0')
    {
      if (cntsl[icnt] == '/')
      {
        hid_t tmp_obj_id;

        // temporarily create a sub-group name,
        cntsl[icnt] = '\0';

        // does the group exist in the file? if not create it
        tmp_obj_id = H5Oopen(rf->h5_file, objname, H5P_DEFAULT);
        if (tmp_obj_id < 0)
          tmp_obj_id = H5Gcreate (rf->h5_file, objname,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // close the group nonetheless
        H5Gclose (tmp_obj_id);

        // restore previous name
        cntsl[icnt] = '/';
      }
      icnt++;
    }
  }

  return;
}


// write method for rlab's atomic datatypes:
// MDR
int
h5_write_rlab_atomic_mdr (Rfile *rf, char *obj_name, MDR *w)
{
  int ires=0;
  hid_t dataset1=-1, fid1=-1;

  if (w->ncol * w->nrow > 1)
  {
    hsize_t fdim[] = {w->ncol, w->nrow};
    fid1 = H5Screate_simple (2, fdim, NULL);
  }
  else if (w->ncol * w->nrow == 1)
    fid1 = H5Screate (H5S_SCALAR);
  else
    return ires;

  // prepare groups in HDF5 file for receipt of dataset
  h5_create_groups_for_writing_atomic_object (rf, obj_name);

  if (w->type == RLAB_TYPE_INT32)
  {
    dataset1 = H5Dcreate (rf->h5_file, obj_name, H5T_NATIVE_INT, fid1,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset1 > -1)
      ires += H5Dwrite(dataset1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, MDIPTR(w));
    else
      ires -= 1;

  }
  else if (w->type == RLAB_TYPE_DOUBLE)
  {
    dataset1 = H5Dcreate (rf->h5_file, obj_name, H5T_NATIVE_DOUBLE, fid1,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    if (dataset1 > -1)
      ires += H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, MDRPTR(w));
    else
      ires -= 1;

  }
  else if (w->type == RLAB_TYPE_FLOAT)
  {
    dataset1 = H5Dcreate (rf->h5_file, obj_name, H5T_NATIVE_FLOAT, fid1,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset1 > -1)
      ires += H5Dwrite(dataset1, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, MDRPTR(w));
    else
      ires -= 1;
  }
  else
    ires = -1;

  // now, add an attribute to the dataset1 if it is not a scalar, which indicates
  // that the dataset has been transposed during writing operation
  if (ires > -1 && (w->ncol)*(w->nrow) > 1)
  {
    int ival = HDF5_RLAB_TRANSPOSE_VAL;
    hid_t space = H5Screate (H5S_SCALAR);
    hid_t attr = H5Acreate  (dataset1, HDF5_RLAB_TRANSPOSE_NAME, H5T_NATIVE_INT, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    ires += H5Awrite (attr, H5T_NATIVE_INT, &ival);
    ires += H5Sclose (space);
    ires += H5Aclose (attr);
  }

  // clean-up
  if (dataset1 > -1)
    H5Dclose (dataset1);

  H5Sclose(fid1);

  return ires;
}

// MDC
int
    h5_write_rlab_atomic_mdc (Rfile *rf, char *obj_name, MDC *wc)
{
  int ires=0;
  hid_t fid1, dataset1, s1_tid;


  // get the data:
  if (wc->ncol * wc->nrow > 1)
  {
    hsize_t fdim[] = {wc->ncol, wc->nrow};
    fid1 = H5Screate_simple (2, fdim, NULL);
  }
  else
    fid1 = H5Screate (H5S_SCALAR);

  // create complex compound
  s1_tid = H5Tcreate (H5T_COMPOUND, sizeof(Complex));
  H5Tinsert(s1_tid, "real", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert(s1_tid, "imag", sizeof(double), H5T_NATIVE_DOUBLE);

  // prepare groups in HDF5 file for receipt of dataset
  h5_create_groups_for_writing_atomic_object (rf, obj_name);

  // write it and close it
  dataset1 = H5Dcreate(rf->h5_file, obj_name, s1_tid, fid1,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset1 > -1)
    ires = H5Dwrite(dataset1, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, MDCPTR(wc));
  else
    ires = -1;

  // now, add an attribute to the dataset1 indicating that the dataset has been
  // transposed:
  if (ires > -1 && (wc->ncol)*(wc->nrow) > 1)
  {
    int ival = HDF5_RLAB_TRANSPOSE_VAL;
    hid_t space = H5Screate (H5S_SCALAR);
    hid_t attr = H5Acreate  (dataset1, HDF5_RLAB_TRANSPOSE_NAME, H5T_NATIVE_INT, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    ires += H5Awrite (attr, H5T_NATIVE_INT, &ival);
    ires += H5Sclose (space);
    ires += H5Aclose (attr);
  }

  // clean-up
  if (dataset1 > -1)
    H5Dclose (dataset1);

  H5Sclose(fid1);
  H5Tclose(s1_tid);

  return ires;
}

//
// MSR -> [nr, nc, nnz, order, ia[nr+1], ja[nnz], d[nnz]]
//
int
    h5_write_rlab_atomic_msr (Rfile *rf, char *obj_name, MSR *w)
{

  int ires = 0;
  hvl_t  wdata[7];

  // scalars:
  wdata[0].len = 1;
  wdata[0].p = (void *) &w->nr;
  wdata[1].len = 1;
  wdata[1].p = (void *) &w->nc;
  wdata[2].len = 1;
  wdata[2].p = (void *) &w->nnz;
  wdata[3].len = 1;
  wdata[3].p = (void *) &w->order;

  // arrays:
  wdata[4].len = (w->nr) + 1;
  wdata[4].p = (void *) w->ia;
  wdata[5].len = (w->nnz);
  wdata[5].p = (void *) w->ja;
  wdata[6].len = (w->nnz);
  wdata[6].p = (void *) w->d;

  hid_t int_array = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t d_array   = H5Tvlen_create (H5T_NATIVE_DOUBLE);

  hid_t memtype  = H5Tcreate (H5T_COMPOUND, 7 * sizeof(hvl_t));
  ires = H5Tinsert (memtype, "nr",    0,               int_array);
  ires = H5Tinsert (memtype, "nc",    1*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "nnz",   2*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "order", 3*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "ia",    4*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "ja",    5*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "d",     6*sizeof(hvl_t), d_array);

  // prepare groups in HDF5 file for receipt of dataset
  h5_create_groups_for_writing_atomic_object (rf, obj_name);

  // get the data:
  hid_t fid1 = H5Screate (H5S_SCALAR);

  // write it and close it
  hid_t dataset1 = H5Dcreate(rf->h5_file, obj_name, memtype, fid1,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset1 > -1)
    ires += H5Dwrite(dataset1, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
  else
    ires = -1;

  // clean-up
  H5Tclose(int_array);
  H5Tclose(d_array);
  if (dataset1 > -1)
    H5Dclose(dataset1);
  H5Sclose(fid1);
  H5Tclose(memtype);

  return ires;
}

//
// MSC -> [nr, nc, nnz, order, ia[nr+1], ja[nnz], c[nnz]]
//
int
    h5_write_rlab_atomic_msc (Rfile *rf, char *obj_name, MSC *w)
{

  int ires = 0;
  hvl_t  wdata[7];

  // scalars:
  wdata[0].len = 1;
  wdata[0].p = (void *) &w->nr;
  wdata[1].len = 1;
  wdata[1].p = (void *) &w->nc;
  wdata[2].len = 1;
  wdata[2].p = (void *) &w->nnz;
  wdata[3].len = 1;
  wdata[3].p = (void *) &w->order;

  // arrays:
  wdata[4].len = (w->nr) + 1;
  wdata[4].p = (void *) w->ia;
  wdata[5].len = (w->nnz);
  wdata[5].p = (void *) w->ja;
  wdata[6].len = 2*(w->nnz);
  wdata[6].p = (void *) w->c;

  hid_t int_array = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t d_array   = H5Tvlen_create (H5T_NATIVE_DOUBLE);

  hid_t memtype  = H5Tcreate (H5T_COMPOUND, 7 * sizeof(hvl_t));
  ires = H5Tinsert (memtype, "nr",    0,               int_array);
  ires = H5Tinsert (memtype, "nc",    1*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "nnz",   2*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "order", 3*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "ia",    4*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "ja",    5*sizeof(hvl_t), int_array);
  ires = H5Tinsert (memtype, "c",     6*sizeof(hvl_t), d_array);

  // prepare groups in HDF5 file for receipt of dataset
  h5_create_groups_for_writing_atomic_object (rf, obj_name);

  // get the data:
  hid_t fid1 = H5Screate (H5S_SCALAR);

  // write it and close it
  hid_t dataset1 = H5Dcreate(rf->h5_file, obj_name, memtype, fid1,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset1 > -1)
    ires += H5Dwrite(dataset1, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
  else
    ires = -1;

  // clean-up
  H5Tclose(int_array);
  H5Tclose(d_array);
  if (dataset1 > -1)
    H5Dclose(dataset1);
  H5Sclose(fid1);
  H5Tclose(memtype);

  return ires;
}


// MDS
int
    h5_write_rlab_atomic_mds (Rfile *rf, char *obj_name, MDS *ws)
{
  int ires=0;

  hid_t space;

  hid_t memtype = H5Tcopy (H5T_C_S1);
  H5Tset_size (memtype, H5T_VARIABLE);

  if (ws->ncol * ws->nrow > 1)
  {
    hsize_t fdim[] = {ws->ncol, ws->nrow};
    space = H5Screate_simple (2, fdim, NULL);
  }
  else
    space = H5Screate (H5S_SCALAR);

  // prepare groups in HDF5 file for receipt of dataset
  h5_create_groups_for_writing_atomic_object (rf, obj_name);

  // write dataset
  hid_t dataset1 = H5Dcreate (rf->h5_file, obj_name, memtype, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset1 > -1)
    ires = H5Dwrite (dataset1, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, MDSPTR(ws));
  else
    ires = -1;

  // now, add an attribute to the dataset1 indicating that the dataset has been
  // transposed:
  if (ires > -1 && (ws->ncol)*(ws->nrow) > 1)
  {
    int ival = HDF5_RLAB_TRANSPOSE_VAL;
    hid_t space = H5Screate (H5S_SCALAR);
    hid_t attr = H5Acreate  (dataset1, HDF5_RLAB_TRANSPOSE_NAME, H5T_NATIVE_INT, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    ires += H5Awrite (attr, H5T_NATIVE_INT, &ival);
    ires += H5Sclose (space);
    ires += H5Aclose (attr);
  }

  // clean-up
  if (dataset1 > -1)
    ires = H5Dclose (dataset1);
  ires = H5Sclose (space);
  ires = H5Tclose (memtype);

  return ires;
}

// read the object assuming that its content is
//    integer, double or string
void *
    get_hobject_atomic_data(hid_t work_obj_id, MDR * coord)
{
  // rlab:
  MDR *w=0;
  MDS *ws=0;

  int transpose = 0;

  int i, j, coord_npoints=0, coord_dim_space=0;
  if (coord)
  {
    coord_npoints = (coord->nrow);
    coord_dim_space = (coord->ncol);
  }

  // HDF:
  hid_t      datatype = H5Dget_type(work_obj_id);
  H5T_class_t t_class = H5Tget_class(datatype);
  hid_t     dataspace = H5Dget_space(work_obj_id);
  int           nrank = H5Sget_simple_extent_ndims(dataspace);
  hid_t attr;

  if (coord && (nrank != coord_dim_space))
  {
    fprintf (stderr, "readm: mismatch in dimensions of dataset and coordinates");

    // clean-up
    H5Sclose(dataspace);
    H5Tclose(datatype);
    return (0);
  }

  hid_t memtype;
  hsize_t       *dims = 0;
  hsize_t   total_size = 1;

  if (nrank > 0)
  {
    dims = (hsize_t *) GC_malloc(nrank * sizeof(hsize_t));
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    total_size = dims[0];
    for (i=1; i<nrank;i++)
      total_size *= dims[i];
  }

  if (nrank > 2)
    fprintf(stderr, "readm: (h5) rlab cannot really handle arrays of rank greater than 2\n");

  switch (t_class)
  {
    case H5T_INTEGER:
      if(!coord_npoints)
      {
        // we are reading an entire dataset
        if (nrank == 0)
          w = mdi_Create(1, 1);
        else if (nrank == 1)
          w = mdi_Create(1, total_size);
        else if (nrank == 2)
          w = mdi_Create(dims[1], dims[0]);
        else if (nrank > 2)
          w = mdi_Create(1, total_size);
        H5Dread(work_obj_id, H5T_NATIVE_INT, H5S_ALL, dataspace,
                H5P_DEFAULT, MDIPTR(w));
      }
      else
      {
        // array 'coord' with
        //  coord_dim_space = coord->ncol; and,
        //  coord_npoints   = coord->nrow;
        // contains the coordinates of the points from dataset that we try to read
        int j1;
        w = mdi_Create(1, coord_npoints);
        hsize_t *cs1 = GC_malloc(nrank * sizeof(hsize_t));
        hsize_t marray = 1;
        hid_t mid1 = H5Screate_simple(1, &marray, NULL);
        for (i=0; i<coord_npoints; i++)
        {
          // read a single point from dataset
          //                 MdiV0(w,i) = ...
          if (coord->type == RLAB_TYPE_INT32)
          {
             for (j1=0;j1<nrank;j1++)
               cs1[j1] = (hsize_t) Mdi0(coord, i, j1) - 1;
          }
          else if (coord->type == RLAB_TYPE_DOUBLE)
          {
            for (j1=0;j1<nrank;j1++)
              cs1[j1] = (hsize_t) Mdr0(coord, i, j1) - 1;
          }

        if (H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, cs1) < 0)
          {
            MdiV0(w, i) = -0;
            continue;
          }
          H5Dread(work_obj_id, H5T_NATIVE_INT, mid1, dataspace,
                  H5P_DEFAULT, &MdiV0(w,i) );
        }
        H5Sclose(mid1);
        GC_free (cs1);
      }

      // check if there is an attribute "Transpose" associated with the dataset
      // if it is, DO NOT TRANSPOSE the dataset.
      if (H5Aexists(work_obj_id, HDF5_RLAB_TRANSPOSE_NAME))
      {
        attr = H5Aopen_name (work_obj_id, HDF5_RLAB_TRANSPOSE_NAME);
        H5Aread (attr, H5T_NATIVE_INT, &transpose);
        H5Aclose (attr);
      }

      break;

    case H5T_FLOAT:
      if(!coord_npoints)
      {
        // we are reading an entire dataset
        if (nrank == 0)
          w = mdr_Create(1, 1);
        else if (nrank == 1)
          w = mdr_Create(1, total_size);
        else if (nrank == 2)
          w = mdr_Create(dims[1], dims[0]);
        else if (nrank > 2)
          w = mdr_Create(1, total_size);
        H5Dread(work_obj_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace,
                H5P_DEFAULT, MDRPTR(w));
      }
      else
      {
        // array 'coord' with
        //  coord_dim_space = coord->ncol; and,
        //  coord_npoints   = coord->nrow;
        // contains the coordinates of the points from dataset that we try to read
        int j1;
        w = mdr_Create(1, coord_npoints);
        hsize_t *cs1 = GC_malloc(nrank * sizeof(hsize_t));
        hsize_t marray = 1;
        hid_t mid1 = H5Screate_simple(1, &marray, NULL);
        for (i=0; i<coord_npoints; i++)
        {
          // read a single point from dataset
          //                 MdiV0(w,i) = ...
          if (coord->type == RLAB_TYPE_INT32)
          {
            for (j1=0;j1<nrank;j1++)
              cs1[j1] = Mdi0(coord, i, j1) - 1;
          }
          else if (coord->type == RLAB_TYPE_DOUBLE)
          {
            for (j1=0;j1<nrank;j1++)
              cs1[j1] = (int) Mdr0(coord, i, j1) - 1;
          }
          if (H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, cs1) < 0)
          {
            MdrV0(w, i) = create_nan();
            continue;
          }
          H5Dread(work_obj_id, H5T_NATIVE_DOUBLE, mid1, dataspace,
                  H5P_DEFAULT, &MdrV0(w,i));
        }
        H5Sclose(mid1);
        GC_free (cs1);
      }

      // check if there is an attribute "Transpose" associated with the dataset
      // if it is, DO NOT TRANSPOSE the dataset.
      if (H5Aexists(work_obj_id, HDF5_RLAB_TRANSPOSE_NAME))
      {
        attr = H5Aopen_name (work_obj_id, HDF5_RLAB_TRANSPOSE_NAME);
        H5Aread (attr, H5T_NATIVE_INT, &transpose);
        H5Aclose (attr);
      }

      break;

    case H5T_STRING:

      memtype = H5Tcopy (H5T_C_S1);
      H5Tset_size (memtype, H5T_VARIABLE);
      char **rdata;

      if(!coord_npoints)
      {
        rdata = (char **) GC_malloc_atomic_ignore_off_page (total_size * sizeof (char *));

        // we are reading an entire dataset
        if (nrank == 0)
          ws = mds_Create(1, 1);
        else if (nrank == 1)
          ws = mds_Create(1, total_size);
        else if (nrank == 2)
          ws = mds_Create(dims[1], dims[0]);
        else if (nrank > 2)
          ws = mds_Create(1, total_size);

        H5Dread(work_obj_id, memtype, H5S_ALL, dataspace, H5P_DEFAULT, rdata);

        // copy data back to rlab
        switch (nrank)
        {
          case 0:
            MdsV0(ws, 0) = cpstr(rdata[0]);
            break;

          case 1:
            for (i=0; i<total_size; i++)
              MdsV0(ws, i) = cpstr(rdata[i]);

          case 2:
            for (i=0; i<dims[1]; i++)
              for (j=0; j<dims[0]; j++)
                Mds0(ws, i, j) = cpstr(rdata[j*dims[1] + i]);

          default:
            for (i=0; i<total_size; i++)
              MdsV0(ws, i) = cpstr(rdata[i]);
        }

        // clean up HDF sh-tuff
        H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata);

      }
      else
      {
        // array 'coord' with
        //  coord_dim_space = coord->ncol; and,
        //  coord_npoints   = coord->nrow;
        // contains the coordinates of the points from dataset that we try to read
        rdata = (char **) GC_malloc_atomic_ignore_off_page (coord_npoints * sizeof (char *));
        int j1;
        ws = mds_Create(1, coord_npoints);
        hsize_t *cs1 = GC_malloc(nrank * sizeof(hsize_t));
        hsize_t marray = 1;
        hid_t mid1 = H5Screate_simple(1, &marray, NULL);
        for (i=0; i<coord_npoints; i++)
        {
          // read a single point from dataset
          //                 MdiV0(w,i) = ...
          if (coord->type == RLAB_TYPE_INT32)
          {
            for (j1=0;j1<nrank;j1++)
              cs1[j1] = Mdi0(coord, i, j1) - 1;
          }
          else if (coord->type == RLAB_TYPE_DOUBLE)
          {
            for (j1=0;j1<nrank;j1++)
              cs1[j1] = (int) Mdr0(coord, i, j1) - 1;
          }
          if (H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, cs1) < 0)
          {
            MdsV0(ws, i) = NULL;
            continue;
          }
          // read
          H5Dread(work_obj_id, memtype, mid1, dataspace,
                  H5P_DEFAULT, &rdata[i]);
          // copy to rlab
          MdsV0(ws, i) = cpstr(rdata[i]);
          // reclaim
          H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata[i]);
        }
        H5Sclose(mid1);
        GC_free (cs1);
      }

      // check if there is an attribute "Transpose" associated with the dataset
      // if it is, DO NOT TRANSPOSE the dataset.
      if (H5Aexists(work_obj_id, HDF5_RLAB_TRANSPOSE_NAME))
      {
        attr = H5Aopen_name (work_obj_id, HDF5_RLAB_TRANSPOSE_NAME);
        H5Aread (attr, H5T_NATIVE_INT, &transpose);
        H5Aclose (attr);
      }

      GC_free (rdata);
      H5Tclose (memtype);
      break;

    default:
      break;
  }

  // clean-up
  H5Sclose (dataspace);
  H5Tclose (datatype);
  if (dims)
    GC_free  (dims);

  if (w)
  {
    if (transpose == 0)
      mdr_Transpose_inplace(w);
    return ((void*)w);
  }
  else if (ws)
  {
    if (transpose == 0)
      mds_Transpose_inplace(ws);
    return ((void *)ws);
  }
  else
    return (NULL);
}

//
// read some special compound types:
//    complex as a pair of two doubles or two floats,
//    where the first one is the real, and the second is the imaginery
//    part of the number.
hid_t
    ascertain_hobject_is_compound_special_complex (hid_t work_obj_id)
{
  hid_t rval=-1;

  // get the type
  hid_t      datatype = H5Dget_type(work_obj_id);
  if (datatype == -1)
    return rval;

  // class should be H5T_COMPOUND
  H5T_class_t t_class = H5Tget_class(datatype);
  if (t_class != H5T_COMPOUND)
  {
    H5Tclose(datatype);
    return rval;
  }

  // should have two members:
  if (H5Tget_nmembers(datatype) != 2)
  {
    H5Tclose(datatype);
    return rval;
  }

  // get their names and check that they are numeric atomic types

  // first member
  // name and type
  char   *mname1 = H5Tget_member_name(datatype, 0);
  hid_t   mtype1 = H5Tget_member_type(datatype, 0);

  if (!mname1)
    goto clean1;
  if (mname1[0] != 'r' && mname1[0] !='R')
    goto clean1;

  // type of the first member should be atomic numeric
  if (    (H5Tget_class(mtype1) != H5T_FLOAT)
      &&  (H5Tget_class(mtype1) != H5T_INTEGER)
     )
    goto clean1;

  // second member
  // name and type
  char   *mname2 = H5Tget_member_name(datatype, 1);
  hid_t   mtype2 = H5Tget_member_type(datatype, 1);

  if (!mname2)
    goto clean2;
  if (mname2[0] != 'i' && mname2[0] !='I')
    goto clean2;
  if (    (H5Tget_class(mtype2) != H5T_FLOAT)
      &&  (H5Tget_class(mtype2) != H5T_INTEGER)
     )
    goto clean2;

  // if we got to this point then it must be a compound with two members,
  // where the first name starts with 'r' and the second with 'i', and the
  // two members are atomic numeric types
  rval = datatype;

clean2:
  // clean-up the second member
  H5Tclose(mtype2);
  if (mname2)
    free (mname2);

clean1:
  // clean-up the first member
  H5Tclose(mtype1);
  if (mname1)
    free (mname1);

  // clean-up top compound if not complex compound
  if (rval == -1)
    H5Tclose(datatype);

  return (rval);
}

//
// this function retrieves the data, but should be called only if the result of
//  ascertain_hobject_is_compound_special_complex
// is greater than 0
//
void *
    get_hobject_compound_special_complex (hid_t work_obj_id, hid_t datatype, MDR * coord)
{
  // rlab:
  MDC *w=0;
  int transpose = 0;

  // integer types
  short int *six=0;
  int        *ix=0;
#if defined (__X86_64__)
  long int  *lix=0;
#endif

  // float types
  float *fx=0;
  double *dx=0;
#if (__SIZEOF_LONG_DOUBLE__ != __SIZEOF_DOUBLE__)
  long double *ldx=0;
#endif

  // buffer for HDF5 reads
  unsigned char *buff=0;

  int i, coord_npoints=0, coord_dim_space=0;
  if (coord)
  {
    coord_npoints = (coord->nrow);
    coord_dim_space = (coord->ncol);
  }

  hid_t     dataspace = H5Dget_space(work_obj_id);

  // first member
  hid_t   mtype1 = H5Tget_member_type(datatype, 0);
  size_t   moff1 = H5Tget_member_offset(datatype, 0);
  size_t   size1 = H5Tget_size(mtype1);

  // second member
  hid_t   mtype2 = H5Tget_member_type(datatype, 1);
  size_t   moff2 = H5Tget_member_offset(datatype, 1);
  size_t   size2 = H5Tget_size(mtype2);

  hsize_t *dims=0;
  hsize_t  total_size;
  int      nrank = H5Sget_simple_extent_ndims(dataspace);
  if (nrank > 0)
  {
    if (coord && (nrank != coord_dim_space))
    {
      fprintf (stderr, "readm: mismatch in dimensions of dataset and coordinates");

    // clean-up
    H5Tclose(mtype1);
    H5Tclose(mtype2);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    return (0);
    }

    // get the dimensions of the dataset and figure out the total size
    dims = GC_malloc(nrank * sizeof(hsize_t));

    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    total_size = dims[0];
    for (i=1; i<nrank;i++) total_size *= dims[i];
  }
  else
    total_size = 1;

  // figure out the size of a compound
  hsize_t msize = H5Tget_size(datatype);

  if (nrank > 2)
    fprintf(stderr, "readm: (h5) rlab cannot really handle arrays of rank greater than 2\n");

  if(!coord_npoints)
    coord_npoints = total_size;

  // we are reading an entire dataset into array 'buff'
  buff = (unsigned char *) GC_malloc_atomic_ignore_off_page (coord_npoints * msize);
  if (!buff)
    rerror ("Terrible internal error: readm h5 cannot initialize the buffer");

  if (coord)
  {
    // we need to duplicate the coordinate array
    // array 'coord' with
    //  coord_dim_space = coord->ncol; and,
    //  coord_npoints   = coord->nrow;
    // contains the coordinates of the points from dataset that we try to read
    int j1;
    int istatus=0;
    hsize_t *cs1 = (hsize_t *) GC_malloc(nrank * sizeof(hsize_t));
    hsize_t marray = 1;
    hid_t mid1 = H5Screate_simple(1, &marray, NULL);

    for (i=0; i<coord_npoints; i++)
    {
      // read a single point from dataset
      //                 MdiV0(w,i) = ...
      if (coord->type == RLAB_TYPE_INT32)
      {
        for (j1=0;j1<nrank;j1++)
          cs1[j1] = (hsize_t) Mdi0(coord, i, j1) - 1;
      }
      else if (coord->type == RLAB_TYPE_DOUBLE)
      {
        for (j1=0;j1<nrank;j1++)
          cs1[j1] = (hsize_t) Mdr0(coord, i, j1) - 1;
      }

      istatus = H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, cs1);
      istatus += H5Dread(work_obj_id, datatype, mid1, dataspace,
              H5P_DEFAULT, &buff[i*msize]);
      if (istatus < 0)
      {
        fprintf (stderr, "readm: error occured during reading for an entry with coordinates\n");
        fprintf (stderr, "readm:   [%i", (int) cs1[0]+1);
        for (j1=1;j1<nrank;j1++)
          fprintf (stderr, ",%i", (int) cs1[j1]+1);
        fprintf (stderr, "].\n");
        fprintf (stderr,
                 "readm: (h5) Please check the entry manually and ignore the result of reading.\n");
      }
    }
    H5Sclose(mid1);
    GC_free (cs1);

    if (istatus == -coord_npoints)
    {
      fprintf(stderr, "readm: (h5) reading of complex dataset from file failed entirely!\n");

      // clean-up
      GC_free (buff);
      if (dims)
        GC_free(dims);
      H5Tclose(mtype1);
      H5Tclose(mtype2);
      H5Sclose(dataspace);
      H5Tclose(datatype);
      return (0);
    }
  }
  else
  {
    // we read an entire dataset
    if( H5Dread(work_obj_id, datatype, H5S_ALL, dataspace,
            H5P_DEFAULT, buff) < 0)
    {
      fprintf(stderr, "readm: (h5) reading of complex dataset from file did not succeed!\n");

      // clean-up
      GC_free (buff);
      if (dims)
        GC_free(dims);
      H5Tclose(mtype1);
      H5Tclose(mtype2);
      H5Sclose(dataspace);
      H5Tclose(datatype);
      return (0);
    }
  }

  // rlab array that will store the result
  if (nrank == 0)
    w = mdc_Create(1, 1);
  else if (nrank == 1)
    w = mdc_Create(1, coord_npoints);
  else if (nrank == 2)
    w = mdc_Create(dims[1], dims[0]);
  else if (nrank > 2)
    w = mdc_Create(1, coord_npoints);

  // now convert buff to data, one coordinate at the time
  // first coordinate:
  if (H5Tget_class(mtype1) == H5T_INTEGER)
  {
    switch (size1)
    {
      case sizeof(int):
        for (i=0; i<coord_npoints; i++)
        {
          ix = (int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *ix;
        }
        break;

      case sizeof(short int):
        for (i=0; i<coord_npoints; i++)
        {
          six = (short int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *six;
        }
        break;

#if defined (__X86_64__)
      case sizeof(long int):
        for (i=0; i<coord_npoints; i++)
        {
          lix = (long int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *lix;
        }
        break;
#endif
    }
  }
  else if (H5Tget_class(mtype1) == H5T_FLOAT )
  {
    switch (size1)
    {
      case sizeof(float):
        for (i=0; i<coord_npoints; i++)
        {
          fx = (float *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *fx;
        }
        break;

      case sizeof(double):
        for (i=0; i<coord_npoints; i++)
        {
          dx = (double *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *dx;
        }
        break;

#if (__SIZEOF_LONG_DOUBLE__ != __SIZEOF_DOUBLE__)
      case sizeof(long double):
        for (i=0; i<coord_npoints; i++)
        {
          ldx = (long double *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *ldx;
        }
#endif
    }
  }
  else
  {
    // clean-up
    GC_free (buff);
    if (dims)
      GC_free(dims);
    H5Tclose(mtype1);
    H5Tclose(mtype2);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    rerror ("Terrible Internal Error: (h5) unknown type!");
  }

  // second coordinate:
  if (H5Tget_class(mtype2) == H5T_INTEGER)
  {
    switch (size2)
    {
      case sizeof(int):
        for (i=0; i<coord_npoints; i++)
        {
          ix = (int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *ix;
        }
        break;

      case sizeof(short int):
        for (i=0; i<coord_npoints; i++)
        {
          six = (short int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *six;
        }
        break;
#if defined (__X86_64__)
      case sizeof(long int):
        for (i=0; i<coord_npoints; i++)
        {
          lix = (long int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *lix;
        }
        break;
#endif
    }
  }
  else if (H5Tget_class(mtype2) == H5T_FLOAT )
  {
    switch (size2)
    {
      case sizeof(float):
        for (i=0; i<coord_npoints; i++)
        {
          fx = (float *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *fx;
        }
        break;

      case sizeof(double):
        for (i=0; i<coord_npoints; i++)
        {
          dx = (double *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *dx;
        }
        break;
#if (__SIZEOF_LONG_DOUBLE__ != __SIZEOF_DOUBLE__)
      case sizeof(long double):
        for (i=0; i<coord_npoints; i++)
        {
          ldx = (long double *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *ldx;
        }
#endif
    }
  }
  else
  {
    // clean-up
    GC_free (buff);
    if (dims)
      GC_free(dims);
    H5Tclose(mtype1);
    H5Tclose(mtype2);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    rerror ("Terrible Internal Error: (h5) unknown type!");
  }

      // check if there is an attribute "Transpose" associated with the dataset
      // if it is, DO NOT TRANSPOSE the dataset.
  if (H5Aexists(work_obj_id, HDF5_RLAB_TRANSPOSE_NAME))
  {
    hid_t attr = H5Aopen_name (work_obj_id, HDF5_RLAB_TRANSPOSE_NAME);
    H5Aread (attr, H5T_NATIVE_INT, &transpose);
    H5Aclose (attr);
  }

  GC_free(buff);

  // clean-up
  H5Tclose(mtype1);
  H5Tclose(mtype2);
  H5Sclose(dataspace);
  H5Tclose(datatype);

  if (dims)
    GC_free(dims);

  if (w)
  {
    if (transpose == 0)
      mdc_NcTranspose_inplace (w);
    return w;
  }
  else
    return (NULL);
}

//
// read some special compound types:
//    matrix sparse real
//    MSR -> [nr, nc, nnz, order, ia[nr+1], ja[nnz], d[nnz]]
hid_t
    ascertain_hobject_is_compound_special_msr (hid_t work_obj_id)
{
  hid_t rval=-1;

  int i=0;

  // get the type
  hid_t datatype = H5Dget_type(work_obj_id);
  if (datatype == -1)
  {
    H5Tclose(datatype);
    return rval;
  }

  // class should be H5T_COMPOUND
  H5T_class_t t_class = H5Tget_class(datatype);
  if (t_class != H5T_COMPOUND)
  {
    H5Tclose(datatype);
    return rval;
  }

  // should have seven members:
  int nmemb = H5Tget_nmembers(datatype);
  if (nmemb != 7)
  {
    H5Tclose(datatype);
    return rval;
  }

  if(   H5Tget_member_index(datatype, "nr")   !=0
     || H5Tget_member_index(datatype, "nc")   !=1
     || H5Tget_member_index(datatype, "nnz")  !=2
     || H5Tget_member_index(datatype, "order")!=3
     || H5Tget_member_index(datatype, "ia")   !=4
     || H5Tget_member_index(datatype, "ja")   !=5
     || H5Tget_member_index(datatype, "d")    !=6
    )
    return (rval);

  int j=0;
  for(i=0; i<7 && !j; i++)
  {
    hid_t mtype = H5Tget_member_type(datatype, i);
    if (H5Tget_class(mtype) != H5T_VLEN)
    { j = 1; }
    H5Tclose(mtype);
  }

  if (j)
    H5Tclose(datatype);
  else
    rval = datatype;

  return (rval);
}

//
// this function retrieves the data, but should be called only if the result of
//  ascertain_hobject_is_compound_special_msr
// is greater than 0
//
void *
    get_hobject_compound_special_msr (hid_t work_obj_id, hid_t datatype)
{
  // rlab:
  MSR *w=0;
  hvl_t rdata[7];

  hid_t     dataspace = H5Dget_space(work_obj_id);
  int           nrank = H5Sget_simple_extent_ndims(dataspace);
  if (nrank != 0)
  {
    fprintf (stderr, "readm: (h5) error in retrieving MSR dataset");

    // clean-up
    H5Sclose(dataspace);
    H5Tclose(datatype);
    //
    H5Tclose(datatype);
    return (0);
  }

  hid_t int_array1 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array2 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array3 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array4 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array5 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array6 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t d_array    = H5Tvlen_create (H5T_NATIVE_DOUBLE);

  hid_t memtype  = H5Tcreate (H5T_COMPOUND, 7 * sizeof(hvl_t));
  H5Tinsert (memtype, "nr",    0,               int_array1);
  H5Tinsert (memtype, "nc",    1*sizeof(hvl_t), int_array2);
  H5Tinsert (memtype, "nnz",   2*sizeof(hvl_t), int_array3);
  H5Tinsert (memtype, "order", 3*sizeof(hvl_t), int_array4);
  H5Tinsert (memtype, "ia",    4*sizeof(hvl_t), int_array5);
  H5Tinsert (memtype, "ja",    5*sizeof(hvl_t), int_array6);
  H5Tinsert (memtype, "d",     6*sizeof(hvl_t), d_array);

  if( H5Dread(work_obj_id, memtype, H5S_ALL, dataspace,
      H5P_DEFAULT, rdata) < 0)
  {
    fprintf(stderr, "readm: (h5) reading of MSR from file did not succeed!\n");

    // clean-up
    H5Tclose(int_array1);
    H5Tclose(int_array2);
    H5Tclose(int_array3);
    H5Tclose(int_array4);
    H5Tclose(int_array5);
    H5Tclose(int_array6);
    H5Tclose(d_array);
    // clean up VL arrays
    H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata);
    //
    H5Sclose(dataspace);
    H5Tclose(memtype);
    //
    H5Tclose(datatype);

    return (0);
  }

//   // convert 'rdata' to MSR
   int *nr, *nc, *nnz;
//    int *iord;
   nr  = (int *) rdata[0].p;
   nc  = (int *) rdata[1].p;
   nnz = (int *) rdata[2].p;
//    iord= (int *) rdata[3].p;
   w = msr_Create( *nr, *nc);
   msr_Setup (w, *nnz);
   memcpy (w->ia, rdata[4].p, ((w->nr)+1)*sizeof(int));
   memcpy (w->ja, rdata[5].p, (w->nnz)*sizeof(int));
   memcpy (w->d,  rdata[6].p, (w->nnz)*sizeof(double));

  // kill the three arrays in rdata that have been reassigned to MSR
  H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata);

  // clean-up
  H5Tclose(int_array1);
  H5Tclose(int_array2);
  H5Tclose(int_array3);
  H5Tclose(int_array4);
  H5Tclose(int_array5);
  H5Tclose(int_array6);
  H5Tclose(d_array);
  //
  H5Sclose(dataspace);
  H5Tclose(memtype);
  //
  H5Tclose(datatype);

  if (w)
    return w;
  else
    return (NULL);
}

//
// read some special compound types:
//    matrix sparse complex
//    MSC -> [nr, nc, nnz, order, ia[nr+1], ja[nnz], c[nnz]]
//  where 'c' is of type 'Complex' struct of two double's, real
//  and imaginery part
hid_t
    ascertain_hobject_is_compound_special_msc (hid_t work_obj_id)
{
  hid_t rval=-1;

  int i, j;

  // get the type
  hid_t datatype = H5Dget_type(work_obj_id);
  if (datatype == -1)
  {
    H5Tclose(datatype);
    return rval;
  }

  // class should be H5T_COMPOUND
  H5T_class_t t_class = H5Tget_class(datatype);
  if (t_class != H5T_COMPOUND)
  {
    H5Tclose(datatype);
    return rval;
  }

  // should have seven members:
  int nmemb = H5Tget_nmembers(datatype);
  if (nmemb != 7)
  {
    H5Tclose(datatype);
    return rval;
  }

  if(   H5Tget_member_index(datatype, "nr")   !=0
     || H5Tget_member_index(datatype, "nc")   !=1
     || H5Tget_member_index(datatype, "nnz")  !=2
     || H5Tget_member_index(datatype, "order")!=3
     || H5Tget_member_index(datatype, "ia")   !=4
     || H5Tget_member_index(datatype, "ja")   !=5
     || H5Tget_member_index(datatype, "c")    !=6
    )
    return (rval);

  j=0;
  for(i=0; i<7 && !j; i++)
  {
    hid_t mtype = H5Tget_member_type(datatype, i);
    if (H5Tget_class(mtype) != H5T_VLEN)
    { j = 1; }
    H5Tclose(mtype);
  }

  if (j)
    H5Tclose(datatype);
  else
    rval = datatype;

  return (rval);
}

//
// this function retrieves the data, but should be called only if the result of
//  ascertain_hobject_is_compound_special_msc
// is greater than 0
//
void *
    get_hobject_compound_special_msc (hid_t work_obj_id, hid_t datatype)
{
  // rlab:
  MSC *w=0;
  hvl_t rdata[7];

  hid_t     dataspace = H5Dget_space(work_obj_id);
  int           nrank = H5Sget_simple_extent_ndims(dataspace);
  if (nrank != 0)
  {
    fprintf (stderr, "readm: (h5) error in retrieving MSR dataset");

    // clean-up
    H5Sclose(dataspace);
    H5Tclose(datatype);
    return (0);
  }

  hid_t int_array1 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array2 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array3 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array4 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array5 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t int_array6 = H5Tvlen_create (H5T_NATIVE_INT);
  hid_t d_array   = H5Tvlen_create (H5T_NATIVE_DOUBLE);

  hid_t memtype  = H5Tcreate (H5T_COMPOUND, 7 * sizeof(hvl_t));
  H5Tinsert (memtype, "nr",    0,               int_array1);
  H5Tinsert (memtype, "nc",    1*sizeof(hvl_t), int_array2);
  H5Tinsert (memtype, "nnz",   2*sizeof(hvl_t), int_array3);
  H5Tinsert (memtype, "order", 3*sizeof(hvl_t), int_array4);
  H5Tinsert (memtype, "ia",    4*sizeof(hvl_t), int_array5);
  H5Tinsert (memtype, "ja",    5*sizeof(hvl_t), int_array6);
  H5Tinsert (memtype, "c",     6*sizeof(hvl_t), d_array);

  if( H5Dread(work_obj_id, memtype, H5S_ALL, dataspace,
      H5P_DEFAULT, rdata) < 0)
  {
    fprintf(stderr, "readm: (h5) reading of MSR from file did not succeed!\n");

      // clean-up
    H5Tclose(int_array1);
    H5Tclose(int_array2);
    H5Tclose(int_array3);
    H5Tclose(int_array4);
    H5Tclose(int_array5);
    H5Tclose(int_array6);
    H5Tclose(d_array);
    // clean up VL arrays
    H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata);
    //
    H5Sclose(dataspace);
    H5Tclose(memtype);
    //
    H5Tclose(datatype);
    return (0);
  }

  // convert 'rdata' to MSR
  int *nr, *nc, *nnz;
//   int *iord;
  nr  = (int *) rdata[0].p;
  nc  = (int *) rdata[1].p;
  nnz = (int *) rdata[2].p;
//   iord= (int *) rdata[3].p;
  w = msc_Create( *nr, *nc);
  msc_Setup (w, *nnz);
  memcpy (w->ia, rdata[4].p, ((w->nr)+1)*sizeof(int));
  memcpy (w->ja, rdata[5].p, (w->nnz)*sizeof(int));
  memcpy (w->c,  rdata[6].p, 2*(w->nnz)*sizeof(double));

  // kill the three arrays in rdata that have been reassigned to MSR
  H5Dvlen_reclaim (memtype, dataspace, H5P_DEFAULT, rdata);

  // clean-up
  H5Tclose(int_array1);
  H5Tclose(int_array2);
  H5Tclose(int_array3);
  H5Tclose(int_array4);
  H5Tclose(int_array5);
  H5Tclose(int_array6);
  H5Tclose(d_array);
  //
  H5Sclose(dataspace);
  H5Tclose(memtype);
  //
  H5Tclose(datatype);

  if (w)
    return w;
  else
    return (NULL);
}


//
// establish compound types:
hid_t
    ascertain_hobject_compound_member_is_atomic (hid_t work_obj_id, MDS *compound_lev,
            H5T_class_t * mem_class, int *ic)
{
  hid_t rval=-1;

  // does compound_lev exist?
  if (!compound_lev)
    return rval;

  // does compound_lev have non-zero size
  int nlev = (compound_lev->nrow) * (compound_lev->ncol);
  if (!nlev)
    return rval;

  // get the object type
  hid_t datatype = H5Dget_type(work_obj_id);

  // class should be H5T_COMPOUND
  H5T_class_t t_class = H5Tget_class(datatype);
  if (t_class != H5T_COMPOUND)
    goto clean0;

  (*ic) = H5Tget_member_index(datatype, MdsV0(compound_lev,0));
  if ((*ic) < 0)
    goto clean0;

  if (nlev == 1)
  {
    hid_t mtype1 = H5Tget_member_type(datatype, (*ic));
    *mem_class = H5Tget_class(mtype1);
    // type of the member should be atomic numeric
    if (    (*mem_class == H5T_FLOAT)
        ||  (*mem_class == H5T_INTEGER)
        ||  (*mem_class == H5T_STRING)
       )
    {
      rval = datatype;
    }

    // clean-up the first member
    H5Tclose(mtype1);
  }
  else
  {
  }

clean0:
  // clean-up top compound if not complex compound
  if (rval == -1)
    H5Tclose(datatype);

  return (rval);
}


//
// this function retrieves the data, but should be called only if the result of
//  ascertain_hobject_compound_member_is_atomic
// is greater than 0
//
#if 0
void *
    get_hobject_compoundL1_member_atomic (Hobject * work_obj, hid_t datatype,
                                          int ic,
                                          MDR *coord)
{
  // rlab:
  MDC *w=0;

  // integer types
  short int *six=0;
  int        *ix=0;
  long int  *lix=0;

  // float types
  float *fx=0;
  double *dx=0;
  long double *ldx=0;

  // buffer for HDF5 reads
  unsigned char *buff=0;

  int i, j1, coord_npoints=0, coord_dim_space=0;
  if (coord)
  {
    coord_npoints = (coord->nrow);
    coord_dim_space = (coord->ncol);
  }

  hid_t     dataspace = H5Dget_space(work_obj_id);

  // locate the member
  hid_t   mtypec = H5Tget_member_type  (datatype, ic);
  size_t   moffc = H5Tget_member_offset(datatype, ic);
  size_t   sizec = H5Tget_size         (mtypec);

  // what is the size of the compound member?

  int           nrank = H5Sget_simple_extent_ndims(dataspace);
  if (coord && (nrank != coord_dim_space))
  {
    fprintf (stderr, "readm: mismatch in dimensions of dataset and coordinates");

    // clean-up
    H5Tclose(mtypec);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    return (0);
  }

  // get the dimensions of the dataset and figure out the total size
  hsize_t       *dims = GC_malloc((nrank+1) * sizeof(hsize_t));

  H5Sget_simple_extent_dims(dataspace, dims, NULL);
  hsize_t   total_size = dims[0];
  for (i=1; i<nrank;i++) total_size *= dims[i];

  // figure out the size of a compound
  hsize_t msizec = H5Tget_size(datatype);

  if (nrank > 2)
    fprintf(stderr, "readm: (h5) rlab cannot really handle arrays of rank greater than 2\n");

  if(!coord_npoints)
    coord_npoints = total_size;

  // we are reading an entire dataset into array 'buff' and from
  // it extracting the requested member
  buff = (unsigned char *) GC_malloc(coord_npoints * msizec);
  if (!buff)
    rerror ("Terrible internal error: readm h5 cannot initialize the buffer");

  if (coord)
  {
    // we need to duplicate the coordinate array
    // array 'coord' with
    //  coord_dim_space = coord->ncol; and,
    //  coord_npoints   = coord->nrow;
    // contains the coordinates of the points from dataset that we try to read
    int j1;
    int istatus;
    hsize_t *cs1 = (hsize_t *) GC_malloc((nrank+1) * sizeof(hsize_t));
    hsize_t marray = 1;
    hid_t mid1 = H5Screate_simple(1, &marray, NULL);

    for (i=0; i<coord_npoints; i++)
    {
      // read a single point from dataset
      //                 MdiV0(w,i) = ...
      if (coord->type == RLAB_TYPE_INT32)
      {
        for (j1=0;j1<nrank;j1++)
          cs1[j1] = (hsize_t) Mdi0(coord, i, j1) - 1;
      }
      else if (coord->type == RLAB_TYPE_DOUBLE)
      {
        for (j1=0;j1<nrank;j1++)
          cs1[j1] = (hsize_t) Mdr0(coord, i, j1) - 1;
      }

      istatus = H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, cs1);
      istatus += H5Dread(work_obj_id, datatype, mid1, dataspace,
                         H5P_DEFAULT, &buff[i*msizec]);
      if (istatus < 0)
      {
        fprintf (stderr, "readm: error occured during reading for an entry with coordinates\n");
        fprintf (stderr, "readm:   [%i", cs1[0]+1);
        for (j1=1;j1<nrank;j1++)
          fprintf (stderr, ",%i", cs1[j1]+1);
        fprintf (stderr, "].\n");
        fprintf (stderr,
                 "readm: (h5) Please check the entry manually and ignore the result of reading.\n");
      }
    }
    H5Sclose(mid1);
    GC_free (cs1);

    if (istatus == -coord_npoints)
    {
      fprintf(stderr, "readm: (h5) reading of complex dataset from file failed entirely!\n");

      // clean-up
      GC_free (buff);
      GC_free(dims);
      H5Tclose(mtype1);
      H5Tclose(mtype2);
      H5Sclose(dataspace);
      H5Tclose(datatype);
      return (0);
    }
  }
  else
  {
    // we read an entire dataset
    if( H5Dread(work_obj_id, datatype, H5S_ALL, dataspace,
        H5P_DEFAULT, buff) < 0)
    {
      fprintf(stderr, "readm: (h5) reading of complex dataset from file did not succeed!\n");

      // clean-up
      GC_free (buff);
      GC_free(dims);
      H5Tclose(mtype1);
      H5Tclose(mtype2);
      H5Sclose(dataspace);
      H5Tclose(datatype);
      return (0);
    }
  }

  // rlab array that will store the result
  if (nrank == 1)
    w = mdc_Create(1, coord_npoints);
  else if (nrank == 2)
    w = mdc_Create(dims[1], dims[0]);
  else if (nrank > 2)
    w = mdc_Create(1, coord_npoints);

  // now convert buff to data, one coordinate at the time
  // first coordinate:
  if (H5Tget_class(mtype1) == H5T_INTEGER)
  {
    switch (size1)
    {
      case sizeof(int):
        for (i=0; i<coord_npoints; i++)
        {
          ix = (int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *ix;
        }
        break;

      case sizeof(short int):
        for (i=0; i<coord_npoints; i++)
        {
          six = (short int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *six;
        }
        break;

      case sizeof(long int):
        for (i=0; i<coord_npoints; i++)
        {
          lix = (long int *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *lix;
        }
        break;
    }
  }
  else if (H5Tget_class(mtype1) == H5T_FLOAT )
  {
    switch (size1)
    {
      case sizeof(float):
        for (i=0; i<coord_npoints; i++)
        {
          fx = (float *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *fx;
        }
        break;

      case sizeof(double):
        for (i=0; i<coord_npoints; i++)
        {
          dx = (double *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *dx;
        }
        break;

      case sizeof(long double):
        for (i=0; i<coord_npoints; i++)
        {
          ldx = (long double *) &buff[i*msize + moff1];
          MdcV0r(w,i) = (double) *ldx;
        }
    }
  }
  else
  {
    // clean-up
    GC_free (buff);
    GC_free(dims);
    H5Tclose(mtype1);
    H5Tclose(mtype2);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    rerror ("Terrible Internal Error: (h5) unknown type!");
  }

  // second coordinate:
  if (H5Tget_class(mtype2) == H5T_INTEGER)
  {
    switch (size2)
    {
      case sizeof(int):
        for (i=0; i<coord_npoints; i++)
        {
          ix = (int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *ix;
        }
        break;

      case sizeof(short int):
        for (i=0; i<coord_npoints; i++)
        {
          six = (short int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *six;
        }
        break;

      case sizeof(long int):
        for (i=0; i<coord_npoints; i++)
        {
          lix = (long int *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *lix;
        }
        break;
    }
  }
  else if (H5Tget_class(mtype2) == H5T_FLOAT )
  {
    switch (size2)
    {
      case sizeof(float):
        for (i=0; i<coord_npoints; i++)
        {
          fx = (float *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *fx;
        }
        break;

      case sizeof(double):
        for (i=0; i<coord_npoints; i++)
        {
          dx = (double *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *dx;
        }
        break;

      case sizeof(long double):
        for (i=0; i<coord_npoints; i++)
        {
          ldx = (long double *) &buff[i*msize + moff2];
          MdcV0i(w,i) = (double) *ldx;
        }
    }
  }
  else
  {
    // clean-up
    GC_free (buff);
    GC_free(dims);
    H5Tclose(mtype1);
    H5Tclose(mtype2);
    H5Sclose(dataspace);
    H5Tclose(datatype);
    rerror ("Terrible Internal Error: (h5) unknown type!");
  }


  GC_free(buff);

  // clean-up
  H5Tclose(mtype1);
  H5Tclose(mtype2);
  H5Sclose(dataspace);
  H5Tclose(datatype);

  GC_free(dims);

  if (w)
  {
    mdc_NcTranspose_inplace (w);
    return w;
  }
  else
    return (NULL);
}

#endif