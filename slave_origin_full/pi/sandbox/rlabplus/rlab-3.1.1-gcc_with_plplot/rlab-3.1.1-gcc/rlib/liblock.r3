//
// basic file locking functionality in rlab
//
static(LOCK_RLAB, LOCK_LIST_ID);
LOCK_RLAB = ".lock.";
if (!exist(LOCK_LIST_ID))
{ LOCK_LIST_ID = <<>>; }

lock = function ( fns )
{
  _this_solver = "lock";

  if (!exist(fns))
  { error(_this_solver + ": Filename required. Cannot continue"); }

  if (class(fns)!= "string")
  { error(_this_solver + ": Filenames have to be string vectors. Cannot continue"); }

  if (isempty(fns))
  { error(_this_solver + ": Filenames have to be string vectors. Cannot continue"); }

  if (any(strlen(fns)==0))
  { error(_this_solver + ": Filenames have to be non-zero string vectors. Cannot continue"); }

  // findout pid of shell
  _spid = text(getpid(), "%.0f");

  // separate filename from directory
  for(i in range(fns))
  {
    // lock file is called:
    //    ../location/.../.lock.filename
    fn = fns[i];
    d  = strsplt(fn, "/");
    n  = length(d);
    d[n] = LOCK_RLAB + d[n];
    lock_file = join(d, "/");

    if (!exist(LOCK_LIST_ID.[lock_file]))
    {
      // file does not exist, so we wait until it is available for write
      LOCK_LIST_ID.[lock_file] = fnctl_lock(lock_file);
    }

    touch(lock_file);

  } // for(i in range(fns))

  return 0;
};

unlock = function ( fns )
{
  _this_solver = "unlock";

  if (!exist(fns))
  { error(_this_solver + ": Filename required. Cannot continue"); }

  if (class(fns)!= "string")
  { error(_this_solver + ": Filenames have to be string vectors. Cannot continue"); }

  if (isempty(fns))
  { error(_this_solver + ": Filenames have to be string vectors. Cannot continue"); }

  if (any(strlen(fns)==0))
  { error(_this_solver + ": Filenames have to be non-zero string vectors. Cannot continue"); }

  // findout pid of shell
  _spid = text(getpid(), "%.0f");

  // separate filename from directory
  rval = ones(fns);
  for(i in range(fns))
  {
    fn = fns[i];
    d  = strsplt(fn, "/");
    n  = length(d);
    d[n] = LOCK_RLAB + d[n];
    lock_file = join(d, "/");

    if (exist(LOCK_LIST_ID.[lock_file]))
    {
      // we proceed to unlock the file
      fnctl_lock(LOCK_LIST_ID.[lock_file]);

      // we no longer own it: remove it from memory and storage
      clear(LOCK_LIST_ID.[lock_file]);
      rm(lock_file);

      rval[i] = 1;

    } else {
      // we did not lock the file, so we cannot unlock it
      rval[i] = 0;
    }
  }
  return rval;
};
