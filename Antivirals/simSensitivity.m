normalvariables;

% kappa=1:0.1:15;
% fprintf('Starting Exp-Exp\n');
% for i=1:length(kappa)
%     fprintf('Current kappa: %.2f\n',kappa(i));
%     expexp(i)=rStar(beta,gamma,sigma,ro,tau,'exp',zeta,'exp',1/kappa(i),ita,pi_k,alpha);
% end
% plot(kappa,expexp,'k');
% hold on;
% 
% kappa=1:15;
% for i=1:length(kappa)
% 	fprintf('Current kappa: %d\n',kappa(i));
% 	rstar(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'exp',zeta,'exp',1/kappa(i),3);
% end
% plot(kappa,rstar,'-or');

kappa=1:0.1:15;
% for i=1:length(kappa)
%     fprintf('Current kappa: %.2f\n',kappa(i));
%     expconst(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'const',kappa(i),ita,pi_k,alpha);
% end
% plot(kappa,expconst,'k');
% hold on;
% 
% for i=1:length(kappa)
% 	fprintf('Current kappa: %.2f\n',kappa(i));
% 	expexp(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'exp',kappa(i),ita,pi_k,alpha);
% end
% plot(kappa,expexp,'k');
% hold on;

for i=1:length(kappa)
	fprintf('Current kappa: %.2f\n',kappa(i));
	constexp(i)=rstarconstexp(beta,gamma,sigma,ro,tau,'const',zeta,'exp',kappa(i),ita,3,alpha);
end
plot(kappa,constexp,'k');
hold on;

% for i=1:length(kappa)
% 	fprintf('Current kappa: %.2f\n',kappa(i));
% 	constconst(i)=rStar2(beta,gamma,sigma,ro,tau,'const',zeta,'const',kappa(i),ita,pi_k,alpha);
% end
% plot(kappa,constconst,'k');
% hold on;

kappa=1:10;
for i=1:length(kappa)
	fprintf('Current kappa: %d\n',kappa(i));
	rstar(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'const',zeta,'exp',kappa(i),3);
end
plot(kappa,rstar,'-or');

for i=1:length(kappa)
    fprintf('Current kappa (2): %d\n',kappa(i));
    rs(i)=simH_constexp(beta,gamma,sigma,ro,tau,'const',zeta,'exp',kappa(i),ita,3,alpha);
end
plot(kappa,rs,'-om');