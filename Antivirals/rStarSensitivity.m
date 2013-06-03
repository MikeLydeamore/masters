beta=2;
gamma=1;
sigma=1;
ro=0.3;
tau=0.3;
onDist='exp';
zeta=1;
%kappa=1:5:300;
kappa=0.1:0.1:5;
offDist='const';
%pi_k=[0 0.1 0.3 0.2 0.3 0.1];
pi_k=[0 0 1];
phi_k=[0 0 0 0 0 0];
alpha=1;
ita=0.0;
N=10000;
gam=0;
expconst=zeros(1,length(kappa));
expexp=zeros(1,length(kappa));
constexp=zeros(1,length(kappa));
constconst=zeros(1,length(kappa));
for i=1:length(kappa)
    fprintf('Current kappa: %.2f\n',kappa(i));
    expconst(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'const',kappa(i),ita,pi_k,alpha);
end


fprintf('Starting Exp-Exp\n');
for i=1:length(kappa)
    fprintf('Current kappa: %.2f\n',kappa(i));
    expexp(i)=rStar2(beta,gamma,sigma,ro,tau,'exp',zeta,'exp',kappa(i),ita,pi_k,alpha);
end
% for i=1:length(kappa)
% 	fprintf('Current kappa: .%.2f\n',kappa(i));
% 	simexpexp(i)=simHousehold(beta,gamma,alpha,sigma,tau,ro,'exp',zeta,'exp',kappa(i),3);
% end

fprintf('Starting Const-Exp\n');
for i=1:length(kappa)
    fprintf('Current kappa: %.2f\n',kappa(i));
    constexp(i)=rStar2(beta,gamma,sigma,ro,tau,'const',zeta,'exp',kappa(i),ita,pi_k,alpha);
end
fprintf('Starting Const-Const\n');
for i=1:length(kappa)
    fprintf('Current kappa: %.2f\n',kappa(i));
    constconst(i)=rStar2(beta,gamma,sigma,ro,tau,'const',zeta,'const',kappa(i),ita,pi_k,alpha);
end
plot(kappa,expconst);
hold on;
plot(kappa,expexp,'r');
%plot(kappa,simexpexp,'-ok');
%legend('Exp-Exp','Simulation');
plot(kappa,constexp,'m');
plot(kappa,constconst,'k');
xlabel('\kappa');
ylabel('R*');
legend('Exp-Const','Exp-Exp','Const-Exp','Const-Const'); 

print -dpng rStarComparison02;
    