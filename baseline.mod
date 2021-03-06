var
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
  vw
  aleph
  nu          // 18 preference shock
  a           // 19 technology shock
  xi          // 20 technology shock
  inot
//welfare measurement 
  welfare     // 21 welfare function
  utility     // 22 utility function
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
  ys          // 21 steady state production
  pi_target   // 22 trend inflation
  omega       // 23 average window length

;
  beta     =0.9937;
  alpha    =0.40; 
  b        =0.6403;
  eta      =1.8896;
  gamma    =3.13;
  theta    =9.30;
  epsilon  =6.00;
  alphap   =0.6115;
  alphaw   =0.5858;
  kappaw   =0.5215;
  kappap   =0.5715;
  phii     =0.8624;
  phipi    =1.2791;   // inflation gap reflection parameter
  phiy     =0.2851;   // output gap reflection parameter
  rhoi     =0.6715;
  rhoa     =0.6535;
  rhonu    =0.3897;
      
//trend inflation
  pi_target=0.00; 
%-----------------------------------------------------------------------------------------%
% model sect
%-----------------------------------------------------------------------------------------%
model;
%-----------------------------------------------------------------------------------------%  
% monetray authority
%-----------------------------------------------------------------------------------------%
  [name='(1)Taylor rule']
   inot=inot(-1)^phii*(is*(pi/pis)^phipi*(y/ys)^phiy)^(1-phii)*xi;
   
   i=inot;
  
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
   
  [name='(10)the wage dispersion']
   vw=(1-alphaw)*(w_star/w)^(-theta)+alphaw*(w(-1)/w*pi(-1)^kappaw/pi)^(-theta)*vw(-1);
  
  [name='(11)the wage dispersion']
   aleph=(1-alphaw)*(w_star/w)^(-theta*(1+eta))+alphaw*(w(-1)/w*pi(-1)^kappaw/pi)^(-theta*(1+eta))*aleph(-1);
 
  [name='(12)inflatin evolution']
   pi^(1-epsilon)=(1-alphap)*pi_star^(1-epsilon)+alphap*pi(-1)^(kappap*(1-epsilon));
 
  [name='(13)optimal price decision']
   (pi_star/pi)^(1+alpha*epsilon/(1-alpha))=epsilon/(epsilon-1)*x1/x2;
 
  [name='(14)auxiliary variable x1']
   x1=lamb*mc*y+alphap*beta*pi(1)^(epsilon/(1-alpha))/pi^(kappap*epsilon/(1-alpha))*x1(1);
 
  [name='(15) auxiliary variable x2']
   x2=lamb*y+alphap*beta*pi(1)^(epsilon-1)*pi^(kappap*(1-epsilon))*x2(1);
 
  [name='(16)wage decision']
   w=(1-alpha)*mc*a*(y/a)^(-alpha/(1-alpha));
 
  [name='(17)resource constraint']
   y=c;
 
  [name='(18)technology shock']
   log(a)=rhoa*log(a(-1))+e_a;
 
  [name='(19)perference shock']
   log(nu)=rhonu*log(nu(-1))+e_nu;
 
  [name='(20)money balance demand']
   m^(-gamma)=lamb*(1-1/i);
  
  [name='(21)Fisher equation']
   r=i/pi(1);
     
//welfare measure
 
  [name='(22)utility function']
   utility=log(c-b*c(-1))-chi*n^(1+eta)*aleph/(1+eta);

  [name='(23)welfare function']
   welfare=utility+beta*welfare(1);
end;
 
%-----------------------------------------------------------------------------------------%
%initial value
%-----------------------------------------------------------------------------------------%
steady_state_model;
  nu=1;
  a=1;
  xi=1;
  n=1.15/3;
  pi=(1+pi_target/100)^(1/4);
  i=pi/beta;
  inot=i;
  r=i/pi;
  pi_star=((pi^(1-epsilon)-alphap*pi^(kappap*(1-epsilon)))/(1-alphap))^(1/(1-epsilon));
  vp=(1-alphap)*(pi/pi_star)^(epsilon/(1-alpha))/(1-alphap*pi^(epsilon/(1-alpha)*(1-kappap)));
  mc=(epsilon-1)/epsilon*(pi_star/pi)^(1+alpha*epsilon/(1-alpha))*(1-beta*alphap*pi^(epsilon/(1-alpha)*(1-kappap)))/(1-beta*alphap*pi^((1-kappap)*(epsilon-1)));
  w=(1-alpha)*mc*(n/vp)^(1-alpha)/n*vp;
  w_star=w*((1-alphaw*pi^((1-kappaw)*(theta-1)))/(1-alphaw))^(1/(1-theta));
  vw=(1-alphaw)*(w_star/w)^(-theta)/(1-alphaw*pi^((1-kappaw)*theta));
  aleph=(1-alphaw)*(w_star/w)^(-theta*(1+eta))/(1-alphaw*pi^((1-kappaw)*theta*(1+eta)));
  chi=(theta-1)/theta*w_star/n^((1-alpha)+eta)*(w_star/w)^(theta*eta)*(1-beta*b)/(1-b)*vp^(1-alpha)*(1-beta*alphaw*pi^((1-kappaw)*theta*(1+eta)))/(1-beta*alphaw*pi^((1-kappaw)*(theta-1)));
  y=a*(n/vp)^(1-alpha);
  c=y;
  lamb=nu*(1-beta*b)/(c-b*c);
  m=(lamb*(1-1/i))^(-1/gamma);
  x1=lamb*mc*y/(1-alphap*beta*pi^((1-kappap)*epsilon/(1-alpha)));
  x2=lamb*y/(1-alphap*beta*pi^((1-kappap)*(epsilon-1)));
  h1=chi*n^(1+eta)*(w/w_star)^(theta*(1+eta))/(1-alphaw*beta*pi^((1-kappaw)*theta*(1+eta)));
  h2=lamb*n*(w/w_star)^theta/(1-alphaw*beta*pi^((theta-1)*(1-kappaw)));
  ys=y;  //steady state production
  is=i;  //steady state interest rate
  pis=pi;//steady state inflation
  //welfare measurement 
  utility=log(c-b*c)-chi*n^(1+eta)*aleph/(1+eta);
  welfare=utility/(1-beta); 
  end;
 %-----------------------------------------------------------------------------------------%
resid;
steady;
check;
model_diagnostics;
%-----------------------------------------------------------------------------------------%
shocks;
  var e_nu ;  stderr 0.0296;
end;
stoch_simul(order=2,irf=40,periods=0,nograph);
% //get welfare value
% Welfare_obs=strmatch('welfare',M_.endo_names,'exact');
% Welfares= oo_.dr.ys(Welfare_obs)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(Welfare_obs));
