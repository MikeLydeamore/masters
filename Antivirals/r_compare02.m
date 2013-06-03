beta=2;
gamma=1;
sigma=1;
ro=0.0;
tau=0.3;
zeta=1;
k=3;
alpha=1;
N=100000;
iter=5;
global testIterations
testIterations=0;
j=zeros(1,iter);
tic;
r=r_fzero(beta,gamma,sigma,ro,tau,zeta,k,alpha);
o=toc;
fprintf('r: %.4f - Completed in %.4fs \n',r,o); 
x=0:0.001:50;
avgInf=zeros(1,length(x));
%plot(x,r*x+log(20),'r');
%drawnow;
for i=1:iter
    testIterations=testIterations+1;
    fprintf('Iteration: %2d | ',i);
    tic;
    [totalInfected, eventTime]=simulateAV(beta,gamma,sigma,alpha,zeta,ro,tau,k,N);
    totalTime=cumsum(eventTime);
    %Work out average
    avgInf(1)=avgInf(1)+20; %start with 20 infected
    for l=1:length(x)
        interval=find((x(l)>=totalTime));
        interval=interval(end);
        %leq=find(totalTime(l)<=x);
        %interval=intersect(geq,leq);
        avgInf(l)=avgInf(l)+totalInfected(interval);
    end
        
    j(i)=toc;
    fprintf('Completed in: %.4fs \n',j(i));
end
fprintf('Total simulation time: %.4fs \n',sum(j));
approxSlope=(log(avgInf(find(x==15))/10)-log(avgInf(find(x==10))/10))/5;
difference=abs(r-approxSlope);
percentDiff=difference/((approxSlope+r)/2)*100;
fprintf('Approximated slope: %.4f \n Error in slope: %.4f (%.2f%%)\n',approxSlope,abs(approxSlope-r),percentDiff);
%plot(x,r*x+log(20),'r','LineWidth',2);
save simAV02.mat;
    