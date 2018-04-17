clear all;close all
load('3s3pGusSSPIFM.mat')
[alpha,beta,v]=butcher2shuosher(A,b,r);
S=[A;b];
c=sum(S,2);
% % Equation 
% u_t = - u_x + nu*u_xx % Here nu is your epsilon
% spacial discreatization
n=50;%201;
x=linspace(0,1,n)'; dx=x(2)-x(1);

% % Build Differentition Matricies
% u_t = -D1u + nu*D2u
%--first-order derivative
D1= diag(ones(n,1))+ diag(-ones(n-1,1),-1);
D1(1,end)=-1;
D1=(1/dx)*D1;

%--second-order derivative
D2= -2*diag(ones(n,1))+ diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
D2(end,1)=1;D2(1,end)=1;
D2=(1/dx)^2*D2;

% %Initial Condition
u0= x>=0.5 & x<1; % Step Function Initial Condition

% Split into linear and 'nonlinear'(or whatever you are not treating exactly)
nu = 1e-1;
L = nu*D2;  % Linear
N = @(u) (-D1)*u; % Nonlinear or not stiff part

Lambda=0.1*[1,.8,.6];

numlamb=length(Lambda);
DT=zeros(1,numlamb);
error=zeros(1,numlamb);

for kk=1:numlamb;
    dt = Lambda(kk)*dx;
    DT(kk)=dt;
    t = 0;Tfinal = 0.5;
    
    un=u0;
    uexact=u0;

    while t < Tfinal;
        dt = min(dt, abs(Tfinal - t));
        t = t + dt ;
        
%--Uncomment method you want--Comment method you do not want---------------   
        % IF Forward Euler
%         un=expm(L*dt)*(un + dt*N(un));
 
        % IFRK 22
%         u1=expm(L*dt)*(un + dt*N(un));
%         un=(1/2)*expm(L*dt)*un + (1/2)*(u1 + dt*N(u1));
        
%         % SSPIFRK33
%         u1=(1/2)*expm((2/3)*L*dt)*un + (1/2)*expm((2/3)*L*dt)*un + (2/3)*expm((2/3)*L*dt)*dt*N(un);
%         u2=(2/3)*expm((2/3)*L*dt)*un + (1/3)*(u1 + (4/3)*dt*N(u1));
%         un=(59/128)*expm(L*dt)*un + (15/128)*expm(L*dt)*(un + (4/3)*dt*N(un)) ...
%             + (27/64)*expm((1/3)*L*dt)*(u2 + (4/3)*dt*N(u2));
        
%    fixing this     % SSPIFRK33
         u1=(v(2))*expm(c(2)*L*dt)*un + (alpha(2,1))*expm(c(2)*L*dt)*un + (beta(2,1))*expm(c(2)*L*dt)*dt*N(un);
         u2=(v(3))*expm(c(3)*L*dt)*un + (alpha(3,2))*expm((c(3)-c(2))*L*dt)*u1 + (beta(3,2))*expm((c(3)-c(2))*L*dt)*dt*N(u1);
         un=(v(4))*expm(L*dt)*un + (alpha(4,1))*expm(L*dt)*un + (beta(4,1))*expm(L*dt)*dt*N(un) ...
             + (alpha(4,3))*expm((c(4)-c(3))*L*dt)*u2 + (beta(4,3))*expm((c(4)-c(3))*L*dt)*dt*N(u2);

%--------------------------------------------------------------------------

         uexact=expm(dt*(-D1 + nu*D2))*uexact; %Tracking exact solution
%         
%         % Plot current solution
%         plot(x,un,'-.')
%         hold on
%         plot(x,uexact,'ro')
%         hold off
%         xlabel('x','fontsize',16); ylabel('u','fontsize',16);
%         title(sprintf('t = %f\n',t));
%         axis([0, 1, 0, 1]);
%         %   %  grid on;
%         pause(.1)
    end

    uexact=expm((-D1+nu*D2)*t)*u0;
    error(kk)=norm(un-uexact,inf);
    
end
Order = diff(log(error))./ diff(log(DT));
order=[nan,Order];

T1 = table(DT',error',order','VariableNames',{'dt' 'Error' 'Order'})

figure(2)
plot(log10(DT),log10(error),'--b','LineWidth',8,'markersize',10)
set(gca,'FontSize',15,'fontweight','b')
xlabel('dt','fontsize',20); ylabel('error','fontsize',20);
title('Order of Convergence','fontsize',20)