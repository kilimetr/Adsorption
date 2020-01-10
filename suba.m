function a = suba(temp)
%
% Return coeficient A as a function of temperature
a1 = 3.646;
a2 = 3.019;
t1 = 4;
t2 = 60; % we use linear interpolation for temperatures between 4 and 60 C

a = a1 + (a2-a1)/(t2-t1) * (temp-t1);
