beta=2;
gamma=1;
alpha=1;
k=4;
N=500;
clf;
% hold off;

% hold on;
% %axis([0 16 0 200])
x=0:0.01:6;
% plot(x,exp(r*x),'r');
hold off;
r=r_fzero(beta,gamma,k,alpha);
[avgInfected, t]=massSimSIR(beta,gamma,alpha,N,k);

subplot(2,1,1), plot(t,avgInfected/100); hold on; plot(x,exp(r*x),'r'); axis([0 30 0 max(avgInfected/100)+1]);

%logInfected
logInf=zeros(1,length(avgInfected));
%Need to ignore zero entires, only change non-zero
change=find(avgInfected>0);
logInf(change)=log(avgInfected(change)/100);

subplot(2,1,2), plot(t,logInf); hold on; plot(x,r*x,'r'); axis([0 30 0 max(logInf)+1]);