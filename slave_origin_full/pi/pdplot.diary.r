// RLaB diary file: pdplot.diary.r. Opened Tue Aug 14 01:29:23 2018

a=[1,2,3;3,1,2]';
writem("|pd -p ../../tmp/junkpd",a);
# it tries to write to home directory and ignores absolute path
system("ps2pdf /tmp/junkpd.ps /tmp/junkpd.pdf");
diary();
