beta=2;
gamma=1;
sigma=1;
ro=0.3;
tau=0.3;
onDist='exp';
zeta=1;
kappa=0.2;
offDist='const';
pi_k=[0 0.1 0.3 0.2 0.3 0.1];
phi_k=[0 0 0 0 0 0];
alpha=1;
ita=0.0;
N=10000;
iter=10;
global testIterations
testIterations=0;
j=zeros(1,iter);
tic;
rs=rStar(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,pi_k,alpha);
p=toc;
tic;
r=r_fzero(beta,gamma,sigma,ro,tau,onDist,zeta,offDist,kappa,ita,pi_k,alpha);
o=toc;
fprintf('r: %.4f - Completed in %.4fs \n',r,o); 
fprintf('R*: %.4f - Completed in %.4fs \n',rs,p);
x=0:0.001:50;
avgInf=zeros(1,length(x));
avgAV=0;
%plot(x,r*x+log(20),'r');
%drawnow;
for i=1:iter
    testIterations=testIterations+1;
    fprintf('Iteration: %2d \n',i);
    tic;
    [totalInfected, eventTime, noAV, numberOfInfected(i), popSize]=sim_expconst(beta,gamma,sigma,alpha,zeta,kappa,ro,tau,ita,pi_k,phi_k,N);
    totalTime=cumsum(eventTime);
    endtime(i)=totalTime(end);
    %Work out average
    avgInf(1)=avgInf(1)+20; %start with 20 infected
    for l=1:length(x)
        interval=find((x(l)>=totalTime));
        interval=interval(end);
        %leq=find(totalTime(l)<=x);
        %interval=intersect(geq,leq);
        avgInf(l)=avgInf(l)+totalInfected(interval);
    end
    avgAV=avgAV+noAV;
        
    j(i)=toc;
    fprintf('Number of AV''s allocated: %d | Epidemic end time: %.4f \n Completed in: %.4fs \n',noAV,endtime(i),j(i));
end
fprintf('Average number of AV''s allocated: %d \n',avgAV/iter);
fprintf('Average number of infections: %d which is %.2f%% of the population \n',sum(numberOfInfected)/10,(sum(numberOfInfected)/10)/popSize * 100);    
fprintf('Total simulation time: %.4fs \n',sum(j));

plot(x,log(avgInf/iter));
axis([0 50 2.5 7.5]);
hold on;
plot(x,r*x+log(12),'r');
    