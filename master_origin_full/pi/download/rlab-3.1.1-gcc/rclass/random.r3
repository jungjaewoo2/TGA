/*
 * rlabplus classes for integer, real, discrete and histogramm random number generators
 *
 *      random.r3
 *
 * copyright 2017, M. Kostrun
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


static(_INIT_);
if (exist(_INIT_))
{
  EOF
}

static(IRNG_MAX_NUM,IRNG_SEED_DEV,IRNG_DEFAULT_NO,IRNG_DEFAULT_NAME);
IRNG_MAX_NUM = 32;
IRNG_DEFAULT_NO = 1;
IRNG_DEFAULT_NAME = "gfsr4";

static(RNG_MAX_NUM,RNG_DEFAULT_NO,RNG_DEFAULT_NAME,RNG_DEFAULT_PARAM);
RNG_DEFAULT_NO = 1;
RNG_DEFAULT_NAME = "uniform";
RNG_DEFAULT_PARAM = [0,1];
RNG_MAX_NUM = 32;

static(DRNG_MAX_NUM,DRNG_DEFAULT_NO);
DRNG_DEFAULT_NO = 1;
DRNG_MAX_NUM = 32;

static(HRNG_MAX_NUM,HRNG_DEFAULT_NO);
HRNG_DEFAULT_NO = 1;
HRNG_MAX_NUM = 32;

hrandom = classdef(arg1,arg2,arg3,arg4)
{
  // static = class-member private declaration:
  //    dynamic memory storage that follows the class-member
  static(_idx, _params_bin, _params_range, _val, _idx_irng);
  _idx  = HRNG_DEFAULT_NO;
  _params_range = [];
  _params_bin = [];
  _idx_irng = IRNG_DEFAULT_NO;

  //
  // process first argument
  //
  if ((type(arg1)=="real") || (type(arg1)=="int"))
  {
    if ((arg1 >= 1) && (arg1 <=HRNG_MAX_NUM))
    {
      _idx = arg1;
    }
  }
  else
  {
    // shift arguments to the right, will use default _idx
    if (type(arg1)!="UNDEF")
    {
      arg4=arg3;
      arg3=arg2;
      arg2=arg1;
    }
  }

  if (type(arg2) == "list")
  {
    if (exist(arg2.bin) && exist(arg2.range))
    {
      if (length(arg2.bin)+1 == length(arg2.range))
      {
        _params_bin = arg2.bin;
        _params_range = arg2.range;
        arg4 = arg3;
      }
    }
  }
  else if (length(arg2) == length(arg3)+1)
  {
    _params_range = arg2;
    if ((type(arg3)=="real") || (type(arg3)=="int"))
    {
      _params_bin = arg3;
    }
    else
    {
      error("vector 'param_bin' must have real positive entries!\n");
    }
  }
  else
  {
    error("vectors 'param_range' and 'param_bin' have to be of the same length!\n");
  }


  if ((type(arg4)=="real") || (type(arg4)=="int"))
  {
    if ((arg4 >= 1) && (arg1 <=RNG_MAX_NUM))
    {
      _idx_irng = arg4;
    }
  }
  else if (type(arg4)=="list")
  {
    if (exist(arg4.idx))
    { _idx_irng = arg4.idx; }
  }

  // this is now default name
  hrng(_idx, _params_range, _params_bin, _idx_irng);

  // public declaration:
  //    after calling the class constructor these are available as
  //    methods to access and modify class-member static variables
  public (idx,val,descr,params,idx_irng);
  descr = "histogram random number generator";
  idx = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=HRNG_MAX_NUM))
      {
        if (x != _idx)
        {
          _idx = x;
          hrng(_idx, _params_range, _params_bin, _idx_irng);
        }
      }
    }
    return _idx;
  };
  idx_irng = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=IRNG_MAX_NUM))
      {
        if (x != _idx_irng)
        {
          _idx_irng = x;
          hrng(_idx, _params_range, _params_bin, _idx_irng);
        }
      }
    }
    return _idx_irng;
  };
  params = function(x,y)
  {
    if ((length(x)==length(y)+1)&&(length(x)>0))
    {
      _params_range = x;
      _params_bin = y;
      hrng(_idx, _params_range, _params_bin, _idx_irng);
    }
    else if (type(x)=="list")
    {
      if (exist(x.bin) && exist(x.range))
      {
        if (length(x.bin)+1 == length(x.range))
        {
          _params_bin = x.bin;
          _params_range = x.range;
          hrng(_idx, _params_range, _params_bin, _idx_irng);
        }
      }
    }
    return <<bin=_params_bin;range=_params_range>>;
  };
  val = function(a1,a2)
  {
    hrng(_idx); // make this drng default generator
    if ((type(a1)=="real") || (type(a1)=="int"))
    {
      if ((type(a2)=="real") || (type(a2)=="int"))
      {
        return hrand(a1,a2);
      }
      return hrand(a1);
    }
    return hrand();
  };
};

drandom = classdef(arg1,arg2,arg3,arg4)
{
  // static = class-member private declaration:
  //    dynamic memory storage that follows the class-member
  static(_idx, _params_val, _params_wgt, _val, _idx_irng);
  _idx  = DRNG_DEFAULT_NO;
  _params_val = [];
  _params_wgt = [];
  _idx_irng = IRNG_DEFAULT_NO;

  //
  // process first argument
  //
  if ((type(arg1)=="real") || (type(arg1)=="int"))
  {
    if ((arg1 >= 1) && (arg1 <=DRNG_MAX_NUM))
    {
      _idx = arg1;
    }
  }
  else
  {
    // shift arguments to the right, will use default _idx
    if (type(arg1)!="UNDEF")
    {
      arg4=arg3;
      arg3=arg2;
      arg2=arg1;
    }
  }

  if (type(arg2) == "list")
  {
    if (exist(arg2.val) && exist(arg2.wgt))
    {
      if (length(arg2.val) == length(arg2.wgt))
      {
        _params_val = arg2.val;
        _params_wgt = arg2.wgt;
        arg4 = arg3;
      }
    }
  }
  else if (length(arg2) == length(arg3))
  {
    _params_val = arg2;
    if ((type(arg3)=="real") || (type(arg3)=="int"))
    {
      _params_wgt = arg3;
    }
    else
    {
      error("vector 'param_wgt' must have real positive entries!\n");
    }
  }
  else
  {
    error("vectors 'param_val' and 'param_wgt' have to be of the same length!\n");
  }


  if ((type(arg4)=="real") || (type(arg4)=="int"))
  {
    if ((arg4 >= 1) && (arg1 <=DRNG_MAX_NUM))
    {
      _idx_irng = arg4;
    }
  }
  else if (type(arg4)=="list")
  {
    if (exist(arg4.idx))
    { _idx_irng = arg4.idx; }
  }

  // this is now default name
  drng(_idx, _params_val, _params_wgt, _idx_irng);

  // public declaration:
  //    after calling the class constructor these are available as
  //    methods to access and modify class-member static variables
  public (idx,val,descr,params,idx_irng);
  descr = "discrete random number generator";
  idx = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=RNG_MAX_NUM))
      {
        if (x != _idx)
        {
          _idx = x;
          drng(_idx, _params_val, _params_wgt, _idx_irng);
        }
      }
    }
    return _idx;
  };
  idx_irng = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=IRNG_MAX_NUM))
      {
        if (x != _idx_irng)
        {
          _idx_irng = x;
          drng(_idx, _params_val, _params_wgt, _idx_irng);
        }
      }
    }
    return _idx_irng;
  };
  params = function(x,y)
  {
    if ((length(x)==length(y))&&(length(x)>0))
    {
      _params_val = x;
      _params_wgt = y;
      drng(_idx, _params_val, _params_wgt, _idx_irng);
    }
    else if (type(x)=="list")
    {
      if (exist(x.val) && exist(x.wgt))
      {
        if (length(x.val) == length(x.wgt))
        {
          _params_val = x.val;
          _params_wgt = x.wgt;
          drng(_idx, _params_val, _params_wgt, _idx_irng);
        }
      }
    }
    return <<val=_params_val;wgt=_params_wgt>>;
  };
  val = function(a1,a2)
  {
    drng(_idx); // make this drng default generator
    if ((type(a1)=="real") || (type(a1)=="int"))
    {
      if ((type(a2)=="real") || (type(a2)=="int"))
      {
        return drand(a1,a2);
      }
      return drand(a1);
    }
    return drand();
  };
};



irandom = classdef(arg1,arg2,arg3)
{
  // static = class-member private declaration:
  //    dynamic memory storage that follows the class-member
  static(_idx, _name, _seed, _val, _write);
  _idx = IRNG_DEFAULT_NO;
  _name = IRNG_DEFAULT_NAME;
  _seed = _iseed();

  //
  // process first argument
  // 
  if ((type(arg1)=="real") || (type(arg1)=="int"))
  {
    if ((arg1 >= 1) && (arg1 <=IRNG_MAX_NUM))
    {
      _idx = arg1;
    }
  }
  else
  {
    // shift arguments to the right, will use default _idx
    if (type(arg1)!="UNDEF")
    {
      arg3=arg2;
      arg2=arg1;
    }
  }

  if (strlen(arg2)> 0)
  {
    _name = arg2;
  }
  if ((type(arg3)=="real") || (type(arg3)=="int"))
  {
    _seed = arg3;
  }

  // this is now default name
  irng(_idx, _name, _seed);

  // public declaration:
  //    after calling the class constructor these are available as
  //    methods to access and modify class-member static variables
  public (idx,seed,name,val, descr,state);
  descr = "integer random number generator";
  idx = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=IRNG_MAX_NUM))
      {
        if (x != _idx)
        {
          irng(_idx, _name, _seed);
        }
      }
    }
    return _idx;
  };
  seed = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if (x != _seed)
      {
        _seed = x;
        irng(_idx, _name, _seed);
      }
    }
    return _seed;
  };
  name = function(x)
  {
    if (strlen(x)>0)
    {
      if (x != _name)
      {
        _name = x;
        irng(_idx, _name, _seed);
      }
    }
    return _name;
  };
  val = function(a1,a2)
  {
    irng(_idx); // make this irng default generator
    if ((type(a1)=="real") || (type(a1)=="int"))
    {
      if ((type(a2)=="real") || (type(a2)=="int"))
      {
        return irand(a1,a2);
      }
      return irand(a1);
    }
    return irand();
  };
  state = function(a1)
  {
    if (type(a1)=="list")
    {
      if ((type(a1.name)=="string")&&(type(a1.state)=="int"))
      {
        irng_state(_idx,a1);
        _name = a1.name;
        _seed = nan();
      }
    }
    return irng_state(_idx);
  };
};

random = classdef(arg1,arg2,arg3,arg4)
{
  // static = class-member private declaration:
  //    dynamic memory storage that follows the class-member
  static(_idx, _name, _params, _val, _idx_irng);
  _idx  = RNG_DEFAULT_NO;
  _name = RNG_DEFAULT_NAME;
  _params = RNG_DEFAULT_PARAM;
  _idx_irng = IRNG_DEFAULT_NO;

  //
  // process first argument
  //
  if ((type(arg1)=="real") || (type(arg1)=="int"))
  {
    if ((arg1 >= 1) && (arg1 <=RNG_MAX_NUM))
    {
      _idx = arg1;
    }
  }
  else
  {
    // shift arguments to the right, will use default _idx
    if (type(arg1)!="UNDEF")
    {
      arg4=arg3;
      arg3=arg2;
      arg2=arg1;
    }
  }

  if (strlen(arg2)> 0)
  {
    _name = arg2;
  }

  if ((type(arg3)=="real") || (type(arg3)=="int"))
  {
    _params = arg3;
  }

  if ((type(arg4)=="real") || (type(arg4)=="int"))
  {
    if ((arg4 >= 1) && (arg4 <=RNG_MAX_NUM))
    {
      _idx_irng = arg4;
    }
  }
  else if (type(arg4)=="list")
  {
    if (exist(arg4.idx))
    { _idx_irng = arg4.idx; }
  }

  // this is now default name
  rng(_idx, _name, _params, _idx_irng);

  // public declaration:
  //    after calling the class constructor these are available as
  //    methods to access and modify class-member static variables
  public (idx,name,val,descr,params,idx_irng,randomize);
  descr = "real random number generator";
  idx = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=RNG_MAX_NUM))
      {
        if (x != _idx)
        {
          _idx = x;
          rng(_idx, _name, _params, _idx_irng);
        }
      }
    }
    return _idx;
  };
  idx_irng = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if ((x>= 1) && (x<=IRNG_MAX_NUM))
      {
        if (x != _idx_irng)
        {
          _idx_irng = x;
          rng(_idx, _name, _params, _idx_irng);
        }
      }
    }
    return _idx_irng;
  };
  name = function(x)
  {
    if (strlen(x)>0)
    {
      if (x != _name)
      {
        _name = x;
        rng(_idx, _name, _params, _idx_irng);
      }
    }
    return _name;
  };
  params = function(x)
  {
    if ((type(x)=="real") || (type(x)=="int"))
    {
      if (any(any(x != _params)))
      {
        _params = x;
        rng(_idx, _name, _params, _idx_irng);
      }
    }
    return _params;
  };
  val = function(a1,a2)
  {
    rng(_idx); // make this irng default generator
    if ((type(a1)=="real") || (type(a1)=="int"))
    {
      if ((type(a2)=="real") || (type(a2)=="int"))
      {
        return rand(a1,a2);
      }
      return rand(a1);
    }
    return rand();
  };
  randomize = function(arg1,arg2,arg3)
  {
    if (!exist(arg1))
    { return 1; }

    // make this irng default generator
    rng(_idx);

    // randomize 'arg1' using rng: see manual for what parameters do
    randomize(arg1,arg2,arg3);

    return 0;
  };
};


_INIT_ = 1;


