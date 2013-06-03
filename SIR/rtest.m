beta=2;
gamma=1;
alpha=1;
%k=4;
%pi_k=[0 0 0 1];
pi_k=[0.2 0.3 0.3 0.2];
N=10000;
clf;
% hold off;

% hold on;
% %axis([0 16 0 200])
x=0:0.01:15;
% plot(x,exp(r*x),'r');
hold off;
r=r_fzero(beta,gamma,alpha,pi_k);
%[avgInfected, t, simsUsed]=massSimSIR(beta,gamma,alpha,N,k);

plot(x,r*x+log(10),'r'); axis([0 30 0 12]); drawnow; hold on;

%plot each realisation
for i=1:20
    tic;
    fprintf('Current iteration: %d',i);
    [totalInfected, eventTime]=simulateSIR(beta,gamma,alpha,N,pi_k);
    logInf=zeros(1,length(totalInfected));
    change=find(totalInfected>0);
    logInf(change)=log(totalInfected(change));
    %subplot(2,1,1), plot(eventTime,totalInfected); hold on;
    plot(cumsum(eventTime),logInf);
    drawnow;
    time=toc;
    fprintf(' | Time to complete: %.4f \n',time);
    %subplot(2,1,2),plot(eventTime,logInf); hold on;
end
%subplot(2,1,1), plot(x,exp(r*x),'r'); axis([0 30 0 5000]);

%subplot(2,1,2), 


% %logInfected
% logInf=zeros(1,length(avgInfected));
% %Need to ignore zero entires, only change non-zero
% change=find(avgInfected>0);
% logInf(change)=log(avgInfected(change)/simsUsed);
% 
% subplot(2,1,2), plot(t,logInf); hold on; plot(x,r*x,'r'); axis([0 30 0 max(logInf)+1]);

plot(x,r*x+log(10),'r','LineWidth',3); axis([0 30 0 12]); drawnow; hold on;