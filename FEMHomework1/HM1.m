% -d(exp(x)du/dx)/dx=-exp(x)(cos(x)-2sin(x)-xcos(x)-xsin(x))
% 0<=x<=1
% u(0)=0; u(1)=cos(1);
% u(x)=xcos(x);

h=[1/4;1/8;1/16;1/32;1/64;1/128];

num=size(h);
for k=1:num(1)
  a=0;
  b=1;
  n=(b-a)/h(k);

  % build the stiff matrix;
  A=sparse(n+1,n+1);
  for j=1:n+1
    if j<=n
      A(j+1,j)=-1/(h(k)^2)*(exp(j*h(k))-exp((j-1)*h(k)));
    endif
    if j>=2
      A(j-1,j)=-1/(h(k)^2)*(exp((j-1)*h(k))-exp((j-2)*h(k)));
    endif
    if j<=n&&j>=2
      A(j,j)=1/(h(k)^2)*(exp(j*h(k))-exp((j-1)*h(k)))+1/(h(k)^2)*(exp((j-1)*h(k))-exp((j-2)*h(k)));
    endif
  endfor
  A(1,:)=0;
  A(1,1)=1;
  A(n+1,:)=0;
  A(n+1,n+1)=1;

  % build the tight hand side vector;
  % maybe it is much better to use gauss-quadrature to compute the
  % integral;
  b=zeros(n+1,1);
  t1=0.7745966692;
  t2=0;
  t3=-t1;
  c1=5/9;
  c2=8/9;
  c3=5/9;
  for i=2:n
    b(i)=gauss_q((i-2)*h(k),(i-1)*h(k))+gauss_q2((i-1)*h(k),i*h(k));
  endfor
  b(1)=0;
  b(n+1)=cos(1);

  % solve the whole equations
  u0=zeros(n+1,1);
  u=MG_2(A,b);
  x=0:h(k):1;
  figure;
  plot(x, u,'b');
  hold on

  y=x.*cos(x);
  plot(x, y,'r');
  hold on 

  % compute the error 
  max_error=max(abs(y'-u));
  disp(['The max error of h=',num2str(h(k)),' is :',num2str(max_error)]);

  %u2=GS(A,u0,b,10000);
  %plot(x,u2,'g');

  %u3=A\b;
  %plot(x,u3,'y');
  
endfor




