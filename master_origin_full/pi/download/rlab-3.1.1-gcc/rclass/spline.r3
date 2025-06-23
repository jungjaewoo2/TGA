/*
 * rlabplus class for splines
 * 
 *      spline.r3
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

bspline = classdef(arg1,arg2,arg3,arg4)
{
  // static = class-member private declaration:
  //    dynamic memory storage that follows the class-member
  static(_data_x, _data_y, _s, _options, _coef, _degree, _knot, _residual, _status);
  // input:
  _data_y = [];
  _data_x = [];
  _s = nan();
  _options = <<>>;
  // output:
  _coef = [];
  _degree = nan();
  _knot = [];
  _residual = [];
  _status = nan();

  // y:
  if (type(arg1)=="real")
  { _data_y = arg1; }
  else if (type(arg1)=="list")
  {
    if (exist(arg1.val)&&exist(arg1.wgt))
    { _data_y = arg1; }
  }

  // x:
  if (type(arg2)=="real")
  { _data_x = arg2; }
  else if (type(arg2)=="list")
  {
    if (exist(arg2.val)&&exist(arg2.wgt))
    { _data_x = arg2; }
  }

  // s:
  if (type(arg3)=="real")
  {
    _s = arg3;
  }

  // options:
  if (type(arg4)=="list")
  {
    _options = arg4;
  }

  if (all(!isnan(_s)) && !isempty(_data_x) && !isempty(_data_y))
  {
    </_coef;_degree;_knot;_residual;_status/>= bsplinefit(_data_y, _data_x, _s, _options);
  }

  public(x,y,smooth,val,knot,options,update);
  x = function(d)
  {
    if (type(d)=="real")
    { _data_x = d; }
    else if (type(d)=="list")
    {
      if (exist(d.val)&&exist(d.wgt))
      { _data_x = d; }
    }
    return _data_x;
  };
  y = function(d)
  {
    if (type(d)=="real")
    { _data_y = d; }
    else if (type(d)=="list")
    {
      if (exist(d.val)&&exist(d.wgt))
      { _data_y = d; }
    }
    return _data_y;
  };
  smooth = function(d)
  {
    if (type(d)=="real")
    {
      _s = d[1];
    }
    return _s;
  };
  knot = function(d)
  {
    if (type(d)=="real")
    {
      if (length(d)>1)
      { _s = d; }
    }
    return _s;
  };
  val = function(x)
  {
    if ((length(coef)==0) || (length(degree)==0) || (length(knot)==0))
    {
      </_coef;_degree;_knot;_residual;_status/>= bsplinefit(_data_y, _data_x, _s, _options);
    }
    if ((length(coef)==0) || (length(degree)==0) || (length(knot)==0))
    {
      return [];
    }
    return bsplineval(x, <<coef=_coef;degree=_degree;knot=_knot>>);
  };
  options = function(z)
  {
    if (type(z)=="list")
    {
      _m = ["degree", "periodic", "maxi", "tol", "convex"];
      for (_m_i in members(z))
      {
        if (all(strindex(_m, _m_i)==0))
        { continue; }
        _options.[_m_i] = z.[_m_i];
      }
    }
    return _options;
  };
  update = function()
  {
    </_coef;_degree;_knot;_residual;_status/>= bsplinefit(_data_y, _data_x, _s, _options);
    return 0;
  };
};


_INIT_ = 1;

