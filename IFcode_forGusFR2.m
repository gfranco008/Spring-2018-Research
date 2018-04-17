clear all;close all

% % Choose method
% filename='3s3pSSPIFM.mat'
% load(filename)
% % or provide a method
A=[0 0;1 0];b=[1/2 1/2]; % RK2
% A=[0 0 0;1 0 0;.25 .25 0]; b=[1/6 1/6 2/3]; % ShuOsher33

S=[A;b];
c=sum(S,2);
s=length(A);

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

% % Initial Condition
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
    
    Nvec = zeros(n,s);

    while t < Tfinal;
        dt = min(dt, abs(Tfinal - t));
        t = t + dt ;
        
%--General solver for any RK method you provide----------------------------
        Nvec(:,1) = N(un);
        
        % Calculate the intermediate stage values
        for i = 2:s+1
            temp = expm(c(i)*L*dt)*un;
            for j = 1:i-1

               temp = temp + dt*expm((c(i)-c(j))*L*dt)*S(i,j)*Nvec(:,j);

            end
            Nvec(:,i) = N(temp);
        end
        % Combine to get u^n+1
        un=temp;
%--------------------------------------------------------------------------

        uexact=expm(dt*(-D1 + nu*D2))*uexact; %Tracking exact solution
        
        % Plot current solution
        plot(x,un,'.')
        hold on
        plot(x,uexact,'ro')
        hold off
        xlabel('x','fontsize',16); ylabel('u','fontsize',16);
        title(sprintf('t = %f\n',t));
        axis([0, 1, 0, 1]);
        %   %  grid on;
        pause(.1)
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