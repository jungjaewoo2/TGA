//
// libsprann.r
//
// some functions motivated by the pattern recognition software that do
// not require c-function coding.
//
//
// by marijan kostrun
// Versions:
//  1. 2000-2001, for neural networks class,
//  2. 2006, as an add-on to sprannlib.

lnorm = function( list1, list2, ntype )
{
  lnorm=[];
  if (!exist(ntype)){ ntype="2";}
  if ( class(list1)!="list" || class(list1)!="list" )
  {
    printf("\n lnorm: one of the arguments is not a list !");
    return [];
  }
  if (size(list1)!=size(list2))
  {
    printf("\n lnorm: list arguments are not of the same size !");
    return [];
  }
  if ( (ntype.nr*ntype.nc!=1)&&(ntype.nr*ntype.nc!=size(list1)) )
  {
    printf("\n lnorm: size of 'ntype' does not match the size of lists !");
    return [];
  }
  for (i in 1:size(list1))
  {
    if ( (list1.[i].nr!=list2.[i].nr) || (list1.[i].nc!=list2.[i].nc) )
    {
      printf("\n lnorm: lists disagree in size of %i-th element !",i);
      return [];
    }

    if ( type(list1.[i]) != type(list2.[i]) )
    {
      printf("\n lnorm: lists disagree in type of %i-th element !",i);
      return [];
    }

    lnorm=[lnorm, norm( list1.[i]-list2.[i], ntype[min(ntype.nr,i)] ) ];
  }
  return lnorm;
};

tfeature = function(feat, clas, opts)
{
  eps0 = 0;
  func = "sigmoid";
  if ( !exist(feat) ){ return nan();}
  if ( !exist(clas) ){ return nan();}
  N = feat.nr * feta.nc;
  C = clas.nr * clas.nc;
  if (class(opts) == "list")
  {
    if (exist(opts.eps))
    {
      if (opts.eps > 0 && opts.eps < 0){ eps0 = eps;}
    }
    if (exist(opts.func))
    {
      if (opts.func == "tanh" || ...
          opts.func == "sigmoid" || ...
          opts.func == "jsigmoid" || ...
          opts.func == "linear"){ func = opts.func;}
    }
  }
  if (func == "tanh"   ) { rent = (-1 + eps) * ones(C,N);}
  if (func == "sigmoid" || ...
      func == "jsigmoid"){ rent = eps * ones(C,N);}
  if (func == "linear" ) { rent =  -2 * ones(C,N);}
  for (i in 1:N)
  {
    for (j in 1:C)
    {
      if (feat[i] == classes[j]){ rent[j;i] = 1 - eps;}
    }
  }
  return rent;
};

//
// sequential validation of a data set:
//  split set in 2 parts, return i-th point as a test set while the
//  1:i-1 data points serve as a training set
//
seqval = function(dataset, i)
{
  //
  // dataset = <<data;feat>>
  //
  if (!exist(dataset))
  {
    printf("\n seqval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (class(dataset)!="list")
  {
    printf("\n seqval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.data))
  {
    printf("\n seqval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.feat))
  {
    printf("\n seqval: dataset = <<data;feat>> is missing!");
    return [];
  }
  //
  // i, a number
  //
  if (!exist(i))
  {
    printf("\n seqval: missing argument 'i'!");
    return [];
  }
  if ((i==1) || (i>data.nr+1))
  {
    printf("\n seqval: Index out of range!");
    return [];
  }

  //
  // training set always exists
  //
  train = <<>>;
  train.data = data[1:i-1;];
  train.feat = feat[1:i-1;];

  //
  // test set may be empty
  //
  test  = <<>>;
  if (i<=data.nr)
  {
    test.data = data[i;];
    test.feat = feat[i;];
  } else {
    test.data = [];
    test.feat = [];
  }
  rval = <<>>;
  rval.test  = test;
  rval.train = train;
  return rval;
};

//
// cross validation of a data set:
//  split set in S parts, return i-th as a test set while the remaining
//  S-1 setparts go to training set
//
crossval = function(dataset, i, s)
{
  //
  // dataset = <<data;feat>>
  //
  if (!exist(dataset))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (class(dataset)!="list")
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.data))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.feat))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }

  //
  // i, a number in range 1:s
  //
  if (!exist(i))
  {
    printf("\n crossval: missing argument 'i'!");
    return [];
  }
  //
  // s, a number greater than unity
  //
  if (!exist(s))
  {
    printf("\n crossval: missing argument 'i'!");
    return [];
  }
  if ((i<1) || (i>s))
  {
    printf("\n crossval: Index out of range!");
    return [];
  }

  test  = <<>>;
  train = <<>>;
  rval  = <<>>;

  nda = dataset.data.nr;
  ncl = int(nda/S);
  if (i == s)
  {
    test.data   = data[ (i-1)*ncl+1:nda; ];
    test.feat   = feat[ (i-1)*ncl+1:nda; ];
    train.data  = data[ 1:(i-1)*ncl; ];
    train.feat  = feat[ 1:(i-1)*ncl; ];
  }
  if (i == 1)
  {
    test.data  = data[ 1:i*ncl; ];
    test.feat  = feat[ 1:i*ncl; ];
    train.data = data[ i*ncl+1:nda; ];
    train.feat = feat[ i*ncl+1:nda; ];
  }
  if ((i > 1) && (i < s))
  {
    test.data  = data[(i-1)*ncl+1:i*ncl;];
    test.feat  = feat[(i-1)*ncl+1:i*ncl;];
    train.data = data[1:(i-1)*ncl,i*ncl+1:nda;];
    train.feat = feat[1:(i-1)*ncl,i*ncl+1:nda;];
  }

  rval.train = train;
  rval.test  = test;
  return rval;
};

//
// classfeat:
//    extract classes from the feature vector or the full
//    dataset=<<data;feat>>.
//
classfeat = function(f)
{
  if (!exist(f))
  {
    printf("\n classfeat: feature or dataset is missing!");
    return [];
  }
  if (class(f)!="num" && class(f)!="string" && class(f)!="list")
  {
    printf("\n classfeat: feature or dataset is missing!");
    return [];
  }


  if (class(f)=="num" || class(f)=="string")
  {
    sf = sort(f).val;
    nf = sf.nr * sf.nc;
  } else {
    if (exist(f.feat))
    {
      sf = sort(f.feat).val;
      nf = sf.nr * sf.nc;
    } else {
      printf("\n classfeat: incorrect feature or dataset!");
      return [];
    }
  }

  j = 1;
  cla[j] = sf[1];
  cnt[j] = 1;
  for (i in 2:nf)
  {
    if ( sf[i] == cla[j] )
    {
      cnt[j] = cnt[j] + 1;
    } else {
      j++;
      cla[j] = sf[i];
      cnt[j] = 1;
    }
  }

  rval = <<>>;
  rval.classes = cla;
  rval.cnt     = cnt;
  return rval;
};

//
// 
//
sepclass = function(dataset)
{
  //
  // dataset = <<data;feat>>
  //
  if (!exist(dataset))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (class(dataset)!="list")
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.data))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }
  if (!exist(dataset.feat))
  {
    printf("\n crossval: dataset = <<data;feat>> is missing!");
    return [];
  }

  z = classfeat(dataset);
  if (isempty(z.classes))
  {
    printf("\n crossval: dataset has no features !");
    return [];
  }

  //
  //
  //
  nclass = z.classes.nr * z.classes.nc;
  sizeM  = dataset.data.nr;
  datalst = <<>> ;
  datacnt = <<>>;
  for (cc in 1:nclass)
  {
    datalst.[ z.classes[cc] ] = zeros(z.cnt, dataset.data.nc);
    datacnt.[ z.classes[cc] ] = 0;
  }
  for (i in 1:sizeM)
  {
    cc = dataset.feat[i];
    datacnt.[ cc ] = datacnt.[ cc ] + 1;
    datalst.[ cc ][ datacnt.[ cc ]; ] = dataset.data[i;] ;
  }
  return datalst;
};

