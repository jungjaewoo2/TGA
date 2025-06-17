psprint = function(FN){
	// invert PLPLOT colors to print on white background 
	// FN -> Postscript filename.ps
	//
	// Usage: psprint("fig.ps")
	// result: saved fig.pdf
	// [azg] Tue Apr  3 00:54:50 PDT 2018

 _plprint ( "/tmp/xxx.ps" , "psc" );

sleep(0.1);
system("sed -e 's@^%%Page: 1 1@{1 exch sub} settransfer\\n%%Page: 1 1 @' </tmp/xxx.ps > "+FN+";" );
system("ps2pdf "+FN+"&& /bin/rm "+FN);
};
