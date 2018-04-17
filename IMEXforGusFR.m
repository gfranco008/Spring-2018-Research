% clear all; close all;
% % Equation 
% u_t = - u_x + nu*u_xx % Here nu is your epsilon
load('3s3pImexK=0.1.mat'); %method
nu = 1e-4;
lambda=.01:.01:2;

n=6;

x=linspace(0,1,n)'; dx=x(2)-x(1);
u0=@(x) x>=0.5 & x<1; % Step Function Initial Condition

% u_t = -D1u + nu*D2u
%--first-order derivative
D1= diag(ones(n,1))+ diag(-ones(n-1,1),-1);
D1(1,end)=-1;
D1=(1/dx)*D1;

%--second-order derivative
D2= -2*diag(ones(n,1))+ diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
D2(end,1)=1;D2(1,end)=1;
D2=(1/dx)^2*D2;

I=eye(n);
tvdFun = @(u) sum([abs(diff(u)); abs((u(1)-u(end)))]);

DT=lambda*dx;

for kk=1:length(lambda);
    
    dt=DT(kk);
    
    P=[A;b(:)'];%Explicit Coeficients
    Q=[Ahat;bhat(:)'];%Implicit Coeficients

    s=length(P)-1; %Number of stages the method is using
    
    t=0;                %Initial Time
    Tfinal=0.5;%25*dt;       %Final Time
    un=u0(x);
    TV=tvdFun(un);
    while t<Tfinal;
        dt = min(dt, abs(Tfinal - t));
        %t = t + dt ;
        u(:,1)=(I-dt*nu*Q(1,1)*D2)\un;
        for i=1:s; %i represents the stage
            utex(:,i)=-D1*u(:,i);     % explicit (saved as it keeps getting reused)
            utimp(:,i)=nu*D2*u(:,i);  % implicit (saved as it keeps getting reused)
            %%% Actual Time Stepping where u(:,i) represents the ith stage
            
            %This builds each stage of the RK method
            u(:,i+1)=un;
            for j=1:i
                u(:,i+1)=u(:,i+1)+dt*(P(i+1,j)*utex(:,j)+Q(i+1,j)*utimp(:,j)); %This is simply a summation operator, i.e u=u+dt*f(u)
            end
            if i~=s;
                u(:,i+1)=(I-dt*nu*Q(i+1,j+1)*D2)\u(:,i+1);                     %Incorporates Diagonal Implicit Element like a BE corrector step.
            end
        end
        un=u(:,i+1);% Update  Un
        t=t+dt;
        % %
%                   plot(x,un)
%                   axis([0,1,-1,2])
%         %
%                   pause(.1)
        
        TV=[TV,tvdFun(un)];
    end
    Violation(kk)=max([diff(TV),10^-16]);
end

%figure
subplot(3,4,10)
plot(lambda,log10(Violation),'g.-','markersize',25)
title(legends(10),'fontsize',16)