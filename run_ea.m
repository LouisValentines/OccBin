%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.


clear

%addpath('toolkit_files\')
global M_ oo_

modnam = 'baseline_ea';
modnamstar = 'baseline_zlb';



constraint = 'i<1-1/beta';
constraint_relax ='i>1-1/beta';

% Pick innovation for IRFs
irfshock =char('e_a');      % label for innovation for IRFs

% Solve nonlinear model
shockscale = 10;
nperiods = 40;
shockssequence = zeros(nperiods,1);
shockssequence(2) =0.0205;
shockssequence = shockssequence * shockscale;



% Solve model, generate model IRFs
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods);


% unpack the IRFs                          
for i=1:M_.endo_nbr
  eval([deblank(M_.endo_names(i,:)),'_uncdifference=zdatalinear(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_difference=zdatapiecewise(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_ss=zdatass(i);']);
end


% get parameter values
for i=1:Mbase_.param_nbr
    eval([Mbase_.param_names(i,:),'=Mbase_.params(i);'])
end

%% Modify to plot IRFs and decision rules
figure('name','技术冲击');

subplot(2,3,1);
plot(100*(i_difference+i_ss-1),'k-','LineWidth',2); hold on;
plot(100*(i_uncdifference+i_ss-1),'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('水平偏离','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='名义利率';
title(titlename,'FontSize',10);

subplot(2,3,2);
plot(100*(pi_difference+pi_ss-1),'k-','LineWidth',2); hold on;
plot(100*(pi_uncdifference+pi_ss-1),'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('水平偏离','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='通货膨胀率';
title(titlename,'FontSize',10);


subplot(2,3,3);
plot(100*(r_difference+r_ss-1),'k-','LineWidth',2); hold on;
plot(100*(r_uncdifference+r_ss-1),'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('水平偏离','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='实际利率';
title(titlename,'FontSize',10);

subplot(2,3,4);
plot(100*y_difference/y_ss,'k-','LineWidth',2); hold on;
plot(100*y_uncdifference/y_ss,'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('百分比','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='产出';
title(titlename,'FontSize',10);

subplot(2,3,5);
plot(100*n_difference/n_ss,'k-','LineWidth',2); hold on;
plot(100*n_uncdifference/n_ss,'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('百分比','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='劳动';
title(titlename,'FontSize',10);

subplot(2,3,6);
plot(100*w_difference/w_ss,'k-','LineWidth',2); hold on;
plot(100*w_uncdifference/w_ss,'r--','LineWidth',2);
xlabel('季度','FontSize',9);
ylabel('百分比','FontSize',9);
grid on;
set(gca, 'GridLineStyle', '--'); 
titlename='工资';
title(titlename,'FontSize',10);

legend('ZLB','Baseline')
axis tight;


