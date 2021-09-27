%-----------------------------------------------------------------------------------------%
% Some notes and suggestions 2021.07.17
% Apple silicon & MacOS Big Sur 11.4 
% Matlab r2020b & dynare 4.5.7 
% dynareOBC need commercial MILP solver, i.e. gurobi optimizer
% Due to dynareOBC updates later than dynare, may cause some errors, so for 
% run stable not try dynare 4.6x.
%-----------------------------------------------------------------------------------------%
var
//sticky economy
  lamb        // 1  Lagrange multiplier
  c           // 2  household consumption
  y           // 3  output
  n           // 4  hours work
  i           // 5  interest rate
  r           // 6  real interest rate
  m           // 7  real money balance
  pi          // 8  inflation 
  pi_star     // 9  target price
  x1          // 10 auxiliary x1
  x2          // 11 auxiliary x2
  mc          // 12 marginal cost
  w           // 13 real wage
  w_star      // 14 target wage
  h1          // 15 auxiliary h1
  h2          // 16 auxiliary h2
  vp          // 17 price dispersion
  nu          // 18 preference shock
  a           // 19 technology shock
  xi          // 20 technology shock
//flexible economy 
  wf          // 21 flexible wage
  yf          // 22 flexible production
  nf          // 23 flexible labour
  lambf       // 24 Lagrange multiplier
  x           // 25 output gap
//welfare measurement 
  welfare     // 26 welfare function
  utility     // 27 utility function
  welfare_f   // 28 welfare function in flexible economy
  utility_f   // 29 utility function in flexible economy
  dwelfare    // 30 welfare gap
  inot
;
%-----------------------------------------------------------------------------------------%
% declare exogenous variables
%-----------------------------------------------------------------------------------------%
varexo
  e_i         // 1  monetary policy shock
  e_a         // 2  technology shock
  e_nu        // 3  preference shock
;
%-----------------------------------------------------------------------------------------%
% declare parameters and steady state
%-----------------------------------------------------------------------------------------%
parameters
  beta        // 1  discount factor
  b           // 2  household consumption habit
  eta         // 3  Frisch elasticity
  alpha       // 4  share of labour 
  gamma       // 5  elasticity of money balance
  theta       // 6  elasticity of labour substitution
  epsilon     // 7  elasticity of substitution
  alphap      // 8  Cavlo price sticky
  alphaw      // 9  Cavlo wage sticky
  chi         // 10 negtive utility of labour input
  kappaw      // 11 wage setting index
  kappap      // 12 price setting index
  rhoi        // 13 AR(1) persistence parameter
  rhoa        // 14 AR(1) persistence parameter
  rhonu       // 15 AR(1) persistence parameter
  phii        // 16 policy rate smoothing parameter
  phipi       // 17 Taylor rule inflation
  phiy        // 18 Taylor rule output
  pis         // 19 steady state inflation
  is          // 20 steady state interest rate
  pi_target   // 21 trend inflation
;
%-----------------------------------------------------------------------------------------%
% model sect
%-----------------------------------------------------------------------------------------%
model;
%-----------------------------------------------------------------------------------------%  
% monetray authority
%-----------------------------------------------------------------------------------------%
  [name='(1)Taylor rule']
   inot=inot(-1)^phii*(is*(pi/pis)^phipi*x^phiy)^(1-phii)*xi;
   i=1;
   log(xi)=rhoi*log(xi(-1))+e_i;  
%-----------------------------------------------------------------------------------------%
//sticky economy

 [name='(2)home Euler equation']
   lamb=nu/(c-b*c(-1))-beta*b*nu(1)/(c(1)-b*c);
 
  [name='(3)home Euler equation']
   lamb=beta*lamb(1)*i/pi(1);
 
  [name='(4)optimal wage decision']
   w_star=theta/(theta-1)*h1/h2;
 
  [name='(5)auxiliary h1']
   h1=chi*(w/w_star)^(theta*(1+eta))*n^(1+eta)+alphaw*beta*pi(1)^(theta*(1+eta))/pi^(kappaw*theta*(1+eta))*(w_star(1)/w_star)^(theta*(1+eta))*h1(1);
 
  [name='(6)auxiliary h2']
   h2=lamb*(w/w_star)^theta*n+alphaw*beta*pi(1)^(theta-1)*pi^((1-theta)*kappaw)*(w_star(1)/w_star)^theta*h2(1);
   
  [name='(7)wage index']
   w^(1-theta)=(1-alphaw)*w_star^(1-theta)+pi^(theta-1)*pi(-1)^(kappaw*(1-theta))*alphaw*w(-1)^(1-theta);
 
  [name='(8)production technology']
   y=a*(n/vp)^(1-alpha);
 
  [name='(9)the price dispersion']
   vp=(pi_star/pi)^(-epsilon/(1-alpha))*(1-alphap)+alphap*pi(-1)^(-kappap*epsilon/(1-alpha))*pi^(epsilon/(1-alpha))*vp(-1);
 
  [name='(10)inflatin evolution']
   pi^(1-epsilon)=(1-alphap)*pi_star^(1-epsilon)+alphap*pi(-1)^(kappap*(1-epsilon));
 
  [name='(11)optimal price decision']
   (pi_star/pi)^(1+alpha*epsilon/(1-alpha))=epsilon/(epsilon-1)*x1/x2;
 
  [name='(12)auxiliary variable x1']
   x1=lamb*mc*y+alphap*beta*pi(1)^(epsilon/(1-alpha))/pi^(kappap*epsilon/(1-alpha))*x1(1);
 
  [name='(13) auxiliary variable x2']
   x2=lamb*y+alphap*beta*pi(1)^(epsilon-1)*pi^(kappap*(1-epsilon))*x2(1);
 
  [name='(14)wage decision']
   w=(1-alpha)*mc*a*(y/a)^(-alpha/(1-alpha));
 
  [name='(15)resource constraint']
   y=c;
 
  [name='(16)technology shock']
   log(a)=rhoa*log(a(-1))+e_a;
 
  [name='(17)perference shock']
   log(nu)=rhonu*log(nu(-1))+e_nu;
 
  [name='(18)money balance demand']
   m^(-gamma)=lamb*(1-1/i);
  
  [name='(19)Fisher equation']
   r=i/pi(1);
   
 // flexible economy  
   
  [name='(20)Lagrange multiplier']
   lambf=nu/(yf-b*yf(-1))-beta*b*nu(1)/(yf(1)-b*yf);
   
  [name='(21)flexible labour']
   wf=theta/(theta-1)*chi*nf^eta/lambf;
   
  [name='(22)flexible wage']
   wf=(epsilon-1)/epsilon*((1-alpha)*a*(y/a)^(-alpha/(1-alpha)));
   
  [name='(23)flexible production']
   yf=a*nf^(1-alpha);
   
  [name='(24)outputgap']
   x=y/yf;
   
//welfare measure
 
  [name='(25)utility function']
   utility=log(c-b*c(-1))-chi*n^(1+eta)/(1+eta);

  [name='(26)welfare function']
   welfare=utility+beta*welfare(1);
   
  [name='(27)utility function in flexible economy']
   utility_f=log(yf-b*yf(-1))-chi/(1+eta)*nf^(1+eta);

  [name='(28)welfare function in flexible economy']
   welfare_f=utility_f+beta*welfare_f(1);

  [name='(29)welfare gap']
   dwelfare=welfare_f-welfare;
end;