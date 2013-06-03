constexp=0;constconst=0;expconst=0;expexp=0; rs1=0; rs2=0;
beta=2;
gamma=1;
sigma=1;
ro=0.3;
tau=0.3;
onDist='exp';
zeta=1;
%kappa=1:5:300;
kappa=1:1:150;
offDist='const';
%pi_k=[0 0.1 0.3 0.2 0.3 0.1];
pi_k=[0 0 1];
phi_k=[0 0 0 0 0 0];
alpha=1;
ita=0.0;
N=10000;
gam=0;

% Plot the 'analytic solution'
kappa=0.5:0.1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %.2f\n',kappa(i));
	constexp(i)=rStar2(beta,gamma,sigma,ro,tau,'const',zeta,'exp',kappa(i),ita,pi_k,alpha);
end
subplot(2,2,1)
plot(kappa,constexp,'k');
hold on;
%Plot old simulation
kappa=1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %d\n',kappa(i));
	rs1(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'const',zeta,'exp',kappa(i),3);
end
plot(kappa,rs1,'-or');
legend('Fast result','Simulation');
title('Const-Exp');


kappa=0.5:0.1:15;
for i=1:length(kappa)
	fprintf('Curent kappa: %.2f\n',kappa(i));
	constconst(i)=rStar2(beta,gamma,sigma,ro,tau,'const',zeta,'const',kappa(i),ita,pi_k,alpha);
end
subplot(2,2,2)
plot(kappa,constconst,'k');
hold on;

%Plot old simulation
kappa=1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %d\n',kappa(i));
	rs1(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'const',zeta,'const',kappa(i),3);
end
plot(kappa,rs1,'-or');
legend('Fast result','Simulation');
title('Const-Const');


% Plot the 'analytic solution'
kappa=0.5:0.1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %.2f\n',kappa(i));
	expconst(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'const',kappa(i),ita,pi_k,alpha);
end
subplot(2,2,3)
plot(kappa,expconst,'k');
hold on;
%Plot old simulation
kappa=1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %d\n',kappa(i));
	rs1(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'exp',zeta,'const',kappa(i),3);
end
plot(kappa,rs1,'-or');
legend('Fast result','Simulation');
title('Exp-Const');


% Plot the 'analytic solution'
kappa=0.5:0.1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %.2f\n',kappa(i));
	expexp(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'exp',kappa(i),ita,pi_k,alpha);
end
subplot(2,2,4)
plot(kappa,expexp,'k');
hold on;
%Plot old simulation
kappa=1:15;
for i=1:length(kappa)
	fprintf('Current kappa: %d\n',kappa(i));
	rs1(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'exp',zeta,'exp',kappa(i),3);
end
plot(kappa,rs1,'-or');
legend('Fast result','Simulation');
title('Exp-Exp');