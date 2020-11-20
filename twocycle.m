% u(0)=u(1)=0;
% so xj=jh, j =  1,...,n-1;
% set fx=-12*x^2+6*x;
% ux = (x-1)x^3

n=100;
A=ones(n-1,1);
a=-1*ones(n-2,1);

A=2*diag(A);
A=A+diag(a,-1)+diag(a,1);
h=1/n;
b0=h:h:1-h;
b=-12*b0.^2+6*b0;
A=1/h^2*A;
x0=zeros(n-1,1);
b=b';

tic;
x=GS(A, x0, b, 100);
plot(b0,x,'r');
hold on
disp(['Computing time of G-S methods after 100 iterations:',num2str(toc)]);

tic;
x=GS(A, x0, b, 1000);
plot(b0,x,'o');
hold on
disp(['Computing time of G-S methods after 1000 iterations:',num2str(toc)]);

tic;
x=GS(A, x0, b, 5000);
plot(b0,x,'c');
hold on
disp(['Computing time of G-S methods after 5000 iterations:',num2str(toc)]);

tic;
x=MG_2(A,b);
disp(['Computing time of MG_twocycle method:',num2str(toc)]);
plot(b0, x,'b');
hold on

u=h:h:1-h;
u=u.^4-u.^3;
plot(b0, u,'y');



