function a=rMat(A,x,known,b)
r=x(1);
x=[x; known];
A=A-diag(r*ones(1,length(A)));
for i=1:length(A)
    a(i)=A(i,:)*x(2:end)+b(i);
end

end