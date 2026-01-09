% This normalizes Upsilon st L=1, and K st q=1
%-------------------------------------------------------------
ssstart    = HANK_SSfile(0,1);
Upsilon    = ssstart(end);
Kss        = ssstart(7);
params(11) = Upsilon;
params(34) = Kss;
%-------------------------------------------------------------