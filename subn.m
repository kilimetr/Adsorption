function n = subn(temp)
%
% Return coefficient "n" as a function of temperature
n1 = 3.277;
n2 = 2.428;
t1 = 4;
t2 = 60;

n = n1 + (n2-n1)/(t2-t1)*(temp-t1); % linear interpolation
