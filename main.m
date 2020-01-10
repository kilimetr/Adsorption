
% main program
% [c] = mol/dm3
% [q] = mol/kg

% set model parameters
epse = 0.434;
epsp = 0.57;
rhop = 1.820; % kg/dm3
rhof = 1.000; % kg/dm3
velo = 10;     % cm/min 
cpp  = 0.25*4.1868*1000; % J/kg.C
cpf  = 1.00*4.1868*1000; % J/kg.C
length = 100; % cm
n = 80;
cin = 0.04; % mol/dm3
tin = 60;   % C

% pack parameters
pars = [epse, epsp, rhop, rhof, velo, cpp, cpf, length, n, cin, tin];

% initial conditions in the column

%ccini = 0.0001; % mol/dm3 (cannot be completely zero)
ccini = 0.4; % mol/dm3 (cannot be completely zero)
ttini = 4;     % C
qqini = suba(ttini)*ccini^(1/subn(ttini));      % mol/kg (just guess, will be calculated)

% pack initial conditions
yg   = zeros(1,3*n);
ypg  = zeros(1,3*n);
yfix = zeros(1,3*n);
ypfix= zeros(1,3*n);

k=1;
for i=1:n
  yg  (k+0) = ccini;
  yfix(k+0) = 1;     % fix "cc"
  yg  (k+1) = qqini;
  yfix(k+1) = 0;     % "qq" will be calculated so qq=fce(cc)
  yg  (k+2) = ttini;
  yfix(k+2) = 1;     % fix "tt"
  k = k+3;
end

% calculate consistent initial conditions
[y0, yp0] = decic(@(t,y,yp) model(t,y,yp,pars), 0, yg, yfix, ypg, ypfix);

% check initial conditions are consistent
model(0,y0,yp0,pars)

% integration in time
tend =1*60;                  % min / end simulation
time =linspace(0,tend,1000);
[tt,yy] = ode15i(@(t,y,yp) model(t,y,yp,pars), time, y0, yp0);
  
% draw graphs 
figure

% inflow and outflow concentrations
subplot(2,2,1)
plot (tt,yy(:,1),tt,yy(:,end-2))
xlabel('time')
ylabel('c')
legend('c_in','c_out')

% inflow and outflow temperatures
subplot(2,2,2)
plot (tt,yy(:,3),tt,yy(:,end-0))
xlabel('time')
ylabel('temperature')
legend('t_in','t_out')

% concentration profile
subplot(2,2,3)
plot (1:n,yy(1,1:3:3*n-2),'g',1:n,yy(end,1:3:3*n-2),'m')
legend('first','last')
hold on
plot (1:n,yy(0.2*end/10,1:3:3*n-2),1:n,yy(0.4*end/10,1:3:3*n-2))
plot (1:n,yy(0.6*end/10,1:3:3*n-2),1:n,yy(1*end/10,1:3:3*n-2))
plot (1:n,yy(1.5*end/10,1:3:3*n-2),1:n,yy(2*end/10,1:3:3*n-2))
plot (1:n,yy(3*end/10,1:3:3*n-2),1:n,yy(4*end/10,1:3:3*n-2))
hold off
xlabel('column profile')
ylabel('c')

% temperature profile
subplot(2,2,4)
plot (1:n,yy(1,3:3:3*n-0),'g',1:n,yy(end,3:3:3*n-0),'m')
legend('first','last')
hold on
plot (1:n,yy(0.2*end/10,3:3:3*n-0),1:n,yy(0.4*end/10,3:3:3*n-0))
plot (1:n,yy(0.6*end/10,3:3:3*n-0),1:n,yy(1*end/10,3:3:3*n-0))
plot (1:n,yy(1.5*end/10,3:3:3*n-0),1:n,yy(2*end/10,3:3:3*n-0))
plot (1:n,yy(3*end/10,3:3:3*n-0),1:n,yy(4*end/10,3:3:3*n-0))
hold off
xlabel('column profile')
ylabel('temperature')

% end of main.m



