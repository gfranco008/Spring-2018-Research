clear all;close all

% Equation 
% u_t = - u_x + eps*u_xx % Here e is your epsilon

% spacial discreatization
n = 50;
x = linspace(0, 1, n)'; 
dx = x(2) - x(1);
epso = 1e-1;
lambda = 0.1;

% time discretization
dt = lambda*dx;
t = 0;
tf = 0.5;
% Build Differentition Matricies

% Forward difference
Dm = diag(ones(n, 1)) + diag(-ones(n-1, 1), -1);
Dm(1, end) = -1;
Dm = (1/dx)*Dm;

% Backward difference
Dp = diag(ones(n, 1)) + diag(-ones(n-1, 1), 1);
Dp(end, 1) = -1;

% Second order
Dc = -(1/dx)*Dm*Dp;


% Initial Condition
u0 = x>=0.5 & x<1; % Step Function Initial Condition
un = u0;
uexact = u0;

L = epso*Dc;  % Linear
N = @(u) (-Dm)*u; % Nonlinear or not stiff part

while t < tf;
        dt = min(dt, abs(tf - t));
        t = t + dt ;
        
        un = expm(-L*t)*(un + dt*N(un));
        uexact = expm(dt*(-Dm + epso*Dc))*uexact; %Tracking exact solution
        
        % Plot current solution
        plot(x,un,'-.')
        hold on
        plot(x,uexact,'ro')
        hold off
        xlabel('x','fontsize',16); 
        ylabel('u','fontsize',16);
        title(sprintf('t = %f\n',t));
        axis([0, 1, 0, 1]);
        %   %  grid on;
        pause(.2)
end
