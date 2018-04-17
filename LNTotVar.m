


%--------------------------------------------------------------------------
%% This code solves IMEX methods using a Linear/ Nonlinear aprroach
% In this case we are solving u_t = - u_x + nu*u_xx % Here nu is your epsilon

% Gustavo Franco Reynoso & Leah Isherwood
% March 02, 2018
%--------------------------------------------------------------------------

function [violation, Lambda] = LNTotVar(name,n,Lambda, L, N)
% % Equation
% u_t = - u_x + nu*u_xx % Here nu is your epsilon
% spacial discreatization
x = linspace(0,1,n)';
dx = x(2)-x(1);

% % Build Differentition Matricies
% u_t = -D1u + nu*D2u

% %Initial Condition
u0 = x>=0.5 & x<1; % Step Function Initial Condition

% Split into linear and 'nonlinear'(or whatever you are not treating exactly

%Lambda = [.1:.05:2] ;

numlamb = length(Lambda);
DT = zeros(1,numlamb);
error = zeros(1,numlamb);

% Information needed for osher 3
% Imports and transforms butcher coeff to osher
%load('3s3pGusSSPIFM.mat') % Speacial method for us 
load('3s3pSSPIFM.mat') % Another methid that zeroes everything
[alpha,beta,v] = butcher2shuosher(A,b,r);
S = [A;b];
c = sum(S,2);

% TV function
TV = @(u) sum(abs(diff(u))) + abs(u(end) - u(1) );

violation = [];
for kk=1:numlamb;
    dt = Lambda(kk)*dx;
    DT(kk)=dt;
    t = 0;
    Tfinal = 0.5;
    totv = [];
    
    un = u0;
    uexact = u0;
    
    while t < Tfinal;
        dt = min(dt, abs(Tfinal - t));
        t = t + dt;
        
        %%--Uncomment method you want--Comment method you do not want---------------
        if lower(name) == "euler"
            % IF Forward Euler
            un = expm(L*dt)*(un + dt*N(un));
            totv = [totv, TV(un)];
        elseif lower(name) == "rk2"
            
            % IFRK 22 %%%%%%%%%% Convert to third order %%%%%%%%%%%%%
            u1 = expm(L*dt)*(un + dt*N(un));
            totv = [totv,TV(u1)];
            un = (1/2)*expm(L*dt)*un + (1/2)*(u1 + dt*N(u1));
            totv = [totv, TV(un)];
            
        else
            
            u1 = (v(2))*expm(c(2)*L*dt)*un + (alpha(2,1))*expm(c(2)*L*dt)*un + (beta(2,1))*expm(c(2)*L*dt)*dt*N(un);
            totv = [totv,TV(u1)];
            u2 = (v(3))*expm(c(3)*L*dt)*un + (alpha(3,2))*expm((c(3)-c(2))*L*dt)*u1 + (beta(3,2))*expm((c(3)-c(2))*L*dt)*dt*N(u1);
            totv = [totv,TV(u2)];
            un = (v(4))*expm(L*dt)*un + (alpha(4,1))*expm(L*dt)*un + (beta(4,1))*expm(L*dt)*dt*N(un) ...
                + (alpha(4,3))*expm((c(4)-c(3))*L*dt)*u2 + (beta(4,3))*expm((c(4)-c(3))*L*dt)*dt*N(u2);
            totv = [totv, TV(un)];
        end
        
        
        % Plot current solution
%                 plot(x,un,'-.')
%                 hold on
%                 plot(x,uexact,'ro')
%                 hold off
%                 xlabel('x','fontsize',16); ylabel('u','fontsize',16);
%                 title(sprintf('lam = %f\n',Lambda(kk)));
%                 axis([0, 1, 0, 1]);
%                 %   %  grid on;
% %                 if numlamb/2 <= kk
% %                     pause(.05)
% %                 else
% %                     pause(.005)
% %                 end
%                  totv = [totv, TV(un)];
        %--------------------------------------------------------------------------
    end
       % Plot current solution
%                 plot(x,un,'-.')
%                 hold on
%                 plot(x,uexact,'ro')
%                 hold off
%                 xlabel('x','fontsize',16); ylabel('u','fontsize',16);
%                 title(sprintf('lam = %f\n',Lambda(kk)));
%                 axis([0, 1, 0, 1]);
    violation(kk) = max(10e-16, max(diff(totv)));
end
% if A=="euler"
%     plot(Lambda,log10(violation), 'o')
% else
%     plot(Lambda,log10(violation), 'x')
% end
% title("Total Variation");
% xlabel('Time Step','fontsize',16);
% legend(A)