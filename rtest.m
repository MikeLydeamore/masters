beta=2;
gamma=1;
alpha=1;
k=6;
N=200;

hold off;
r=r_fzero(beta,gamma,k,alpha);
massSimSIR(beta,gamma,alpha,N,k);
hold on;
axis([0 16 0 200])
x=0:0.01:6;
plot(x,exp(r*(x+1.7)),'r');