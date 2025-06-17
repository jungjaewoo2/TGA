//
// this demonstrate rlab <-> firefox connection through Mozlab/repl interface
//
py.init();

host = "localhost";
port = 4242;
pmpt = "repl>"; // this is only valid if python is the only one talking to repl. RTFM
tout = 1;

if (!exist(telnet))
{
	telnet=<<>>;
	telnet.read = function()
	{
		global(pmpt,py,tout);
		_pycmd = "tn.read_until(\"" + pmpt + "\"," + text(tout, "%g") + ")";
		rval = py.eval(_pycmd);
		_i = strindex(rval,pmpt);
		if(_i > 1)
		{ rval = substr(rval, 1:(_i-1)); }

		// use '\n' to separate rval into an array
		_c  = ascii(rval);
		_js = find(_c == 10);
		if (_js.n == 0)
		{
			return rval;
		else
			return strsplt(rval, "\n");
		}
	};

	telnet.write = function( s )
	{
		global(pmpt,py,tout);
		for (_s in s)
		{
			_pycmd = "tn.write(\"" + _s + "\")";
			rval = py.eval(_pycmd);
		}
		return rval;
	};
}

if(!exist(_pycmds1))
{
	_pycmds1 = [ ...
		"import sys", ...
		"from telnetlib import Telnet", ...
		"tn = Telnet('" + host + "'," + text(port, "%g") + ")" ...
	];
	py.cmd( _pycmds1 );

	// read welcome text from 'repl'
	telnet.read();
}

//
// executing some random repl code:
//
telnet.write("document.title");
x = telnet.read();
printf("The window title on your firefox session is:\n%s\n", x);

telnet.write("content.location.href");
x = telnet.read();
printf("Your web session is on the site:\n%s\n", x);

