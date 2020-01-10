function equval = model(time,yvec,ypvec,pars)
% 
% Model for adsorption column
% - no mass or temperature diffusion diffusion
% - instantenous adsorption/desorption
% - no mass/heat transfer resistance between fixed and mobile phase

% extract model parameters
epse = pars(1);     % interparticle porosity (external)
epsp = pars(2);     % intraparticle porosity (particle)
rhop = pars(3);     % density of solid adsorbent
rhof = pars(4);     % density of fluid
velo = pars(5);     % superficial velocity in the column
cpp  = pars(6);     % heat capacity of solid adsorbent
cpf  = pars(7);     % heat capacity of fluid
length = pars(8);   % length of the column
n    = pars(9);     % number of sections (discretization)
cin  = pars(10);    % inflow concentration
tin  = pars(11);    % inflow temperature

dL = length/n; % length of single section

% extract state variables and theirs derivatives
cc = zeros(n,1);
qq = zeros(n,1);
tt = zeros(n,1);
dccdt = zeros(n,1);
dqqdt = zeros(n,1);
dttdt = zeros(n,1);

k=1;
for i=1:n
  cc(i)    = yvec (k+0);
  dccdt(i) = ypvec(k+0);
  qq(i)    = yvec (k+1);
  dqqdt(i) = ypvec(k+1);
  tt(i)    = yvec (k+2);
  dttdt(i) = ypvec(k+2);
  k = k+3;
end 

% calculate residuals
equval = zeros(3*n,1);

% first section
% - mass balance
equval(1) = (epse+(1-epse)*epsp)*dccdt(1) ...
          + (1-epse)*(1-epsp)*rhop*dqqdt(1) ...
          - velo/dL * (cin-cc(1));

% - adsorption equilibrium
ai = suba(tt(1));  % A and N depend on the temperature
ni = subn(tt(1));
equval(2) = qq(1) - ai*cc(1)^(1/ni);

% - enthalpy balance
equval(3) = (epse+epsp*(1-epse))*rhof*cpf*dttdt(1) ...
          + (1-epse)*(1-epsp)*rhop*cpp*dttdt(1) ...
          - velo/dL * cpf*rhof * (tin-tt(1)); 

% now equations for the remaining sections
k=4;
for i=2:n

% - mass balance
  equval(k+0) = (epse+(1-epse)*epsp)*dccdt(i) ...
              + (1-epse)*(1-epsp)*rhop*dqqdt(i) ...
              - velo/dL * (cc(i-1)-cc(i));

% - adsorption equilibrium
  ai = suba(tt(i));
  ni = subn(tt(i));
  equval(k+1) = qq(i) - ai*cc(i)^(1/ni);

% - enthalpy balance
  equval(k+2) = (epse+epsp*(1-epse))*rhof*cpf*dttdt(i) ...
              + (1-epse)*(1-epsp)*rhop*cpp*dttdt(i) ...
              - velo/dL * cpf*rhof * (tt(i-1)-tt(i)); 
  k = k+3;
end

% end of model.m








