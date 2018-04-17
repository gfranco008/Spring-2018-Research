clear all;close all

% % Equation
% u_t = - u_x + nu*u_xx % Here nu is your epsilon
% spacial discreatization
n=51;
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
nu = 1e-2;
L = nu*D2;  % Linear
N = @(u) (-D1)*u; % Nonlinear or not stiff part

Lambda = [.01:.05:2] ;

numlamb=length(Lambda);
DT=zeros(1,numlamb);
error=zeros(1,numlamb);

TV = @(u) sum(abs(diff(u))) + abs(u(1) - u(end));

violation = [];
for kk=1:numlamb;
    dt = Lambda(kk)*dx;
    DT(kk)=dt;
    t = 0;Tfinal = 0.5;
    totv = [];
    
    un=u0;
    uexact=u0;
    
    while t < Tfinal;
        dt = min(dt, abs(Tfinal - t));
        t = t + dt ;
        
        %--Uncomment method you want--Comment method you do not want---------------
        % IF Forward Euler
        un=expm(L*dt)*(un + dt*N(un));
        
        % IFRK 22
        %u1=expm(L*dt)*(un + dt*N(un));
        %un=(1/2)*expm(L*dt)*un + (1/2)*(u1 + dt*N(u1));
        %--------------------------------------------------------------------------
        
        uexact=expm(dt*(-D1 + nu*D2))*uexact; %Tracking exact solution
        
        % Plot current solution
                plot(x,un,'-.')
                hold on
                plot(x,uexact,'ro')
                hold off
                xlabel('x','fontsize',16); ylabel('u','fontsize',16);
                title(sprintf('t = %f\n',t));
                axis([0, 1, 0, 1]);
                %   %  grid on;
                pause(.005)
                 totv = [totv, TV(un)];
    end
    violation(kk) = max(10e-16, max(diff(totv)));
    uexact=expm((-D1+nu*D2)*t)*u0;
    error(kk)=norm(un-uexact,inf);
    
end
% Order = diff(log(error))./ diff(log(DT));
% order=[nan,Order];
%
% T1 = table(DT',error',order','VariableNames',{'dt' 'Error' 'Order'})
%
% figure(2)
% plot(log10(DT),log10(error),'--b','LineWidth',8,'markersize',10)
% set(gca,'FontSize',15,'fontweight','b')
% xlabel('dt','fontsize',20); ylabel('error','fontsize',20);
% title('Order of Convergence','fontsize',20)