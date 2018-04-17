%--------------------------------------------------------------------------
%% This code determines the best distribution for solving IMEX methods using a Linear/ Nonlinear aprroach
% In this case we are solving u_t = - u_x + nu*u_xx % Here nu is your epsilon

% Gustavo Franco Reynoso
% March 02, 2018
%--------------------------------------------------------------------------

clear all;
close all;

n = 6;
x = linspace(0,1,n)';
dx = x(2)-x(1);
nu = 1e-4;

D1 = diag(ones(n,1))+ diag(-ones(n-1,1),-1);
D1(1,end) = -1;
D1 = -(1/dx)*D1;

% % -- 3 points 
% D1 = 2*diag(ones(n,1))+ diag(-(3/2)*ones(n-1,1),-1) + diag(-(1/2)*ones(n-1,1),1);
% D1(1, end) = -3/2;
% D1(end, 1) = -1/2;
% D1 = (1/dx)*D1;


% %--second-order derivative
D2 = -2*diag(ones(n,1))+ diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
D2(end,1) = 1;
D2(1,end) = 1;
D2 = nu*(1/dx)^2*D2;

% % - -- - Centered diff with 4 pnts
% 
% D2 = -(5/2)*diag(ones(n,1)) + (4/3)*diag(ones(n-1,1),1) + (4/3)*diag(ones(n-1,1),-1) + -(1/12)*diag(ones(n-2,1),-2) + -(1/12)*diag(ones(n-2,1),2);
% D2(end,1) = (4/3);
% D2(1,end) = (4/3);
% D2(1,end-1) = -(1/12);
% D2(2,end) = -(1/12);
% D2(end-1,1) = -(1/12);
% D2(end,2) = -(1/12);
% D2 = nu*(1/dx)^2*D2;

% % Split into linear and 'nonlinear'(or whatever you are not treating exactly)
% Divide each matrix into only one diagonal
A = diag(-ones(n-1,1),-1);
A(1,end) = -1;
A = -(1/dx)*A;

B = diag(ones(n,1));
B = -(1/dx)*B;

C = diag(ones(n-1,1),-1);
C(1,end) = 1;
C = nu*(1/dx^2)*C;

D = -2*diag(ones(n,1));
D = nu*(1/dx^2)*D;

E = diag(ones(n-1,1),1);
E(end,1) = 1;
E = nu*(1/dx^2)*E;

fbp = {A, B, C, D, E};

% Save for later
step = [];
steprk = [];
forplot = {};
forplotrk = {};
m = length(fbp);
k = 0;
lam = {};
lamrk = {};

tic()

for i = 1:length(fbp)-1
    comb = nchoosek(1:m,i); % Create all posible combinations choosing 1:m
    j = 0;
    while j < length(comb)
        j = j+1;
        k = k+1;
        comb(j,:); % choose a row
        fbp = {A, B, C, D, E};
        ofbp = fbp(comb(j,:));
        linmat = sum(cat(3,ofbp{:}),3); % sum the chosen matrixes for L
        fbp(comb(j,:)) = [];
        nonmat = sum(cat(3,fbp{:}),3); % sum the chosen matrixes for N
        [violation, Lambda] = LNTotVar("euler", n, linspace(.01, 2), linmat, @(u) nonmat*u);
         
        % This determines the max number of points before there
                    % is a jump in TV
        pnt = min(find(diff(violation)> 1e-4))-1;
        
        forplot{k} = violation;
        lam{k} = Lambda;
        if isempty(pnt)==1
            pnt = 0;
        end
        step = [step; pnt, string(char(comb(j,:) + 64))];
        % hold on
        [violationrk, Lambdark] = LNTotVar("rk3", n, linspace(.01, 2), linmat, @(u) nonmat*u);
        pnt = min(find(diff(violationrk)> 1e-4))-1;
        forplotrk{k} =  violationrk;
        lamrk{k} = Lambdark;
        if isempty(pnt)==1
            pnt = 0;
        end
        
         % Print out the inputs
        steprk = [steprk; pnt, string(char(comb(j,:) + 64))];
        % legend("euler", "rk2")
        % hold off
        % figure()
        %         if j>5
        %             keyboard
        %         end
    end
    
end
indx = find(str2double(step(:,1)) >= max(str2double(step(:,1))));
indrk = find(str2double(steprk(:,1)) >= max(str2double(steprk(:,1))));
pp = [ 's' 'o' '+' '*' 'd' '<' 'h' '-' 's' 'o' '+' '*' 'd' '<' 'h' '-' 's' 'o' '+' '*' 'd' '<' 'h' '-' 's' 'o' '+' '*' 'd' '<' 'h' '-' 's' 'o' '+' '*' 'd' '<' 'h' '-' 's' 'o' '+' '*' 'd' '<' 'h' '-'];
% Chose the best seperation using Euler 
for j = 1:length(indx)
    
    % Convert to matrix
    i = indx(j)
    violation = cell2mat(forplot(i));
    violationrk = cell2mat(forplotrk(i));
    
    
    plot(cell2mat(lam(i)),log10(violation), pp(j))
    %legend(step(indrk,2))
    hold on
    plot(cell2mat(lamrk(i)),log10(violationrk), pp(j))
    
    title("Best Total Variation for Euler",'fontsize',16);
    xlabel(strjoin(["Time Step : ", step(i,2)]), 'fontsize',16);
    
    legend("euler", "rk2")
    
%     savefig(j, char(strjoin(["euler_", step(i,2), ".fig"], "")));
%     saveas(j,char(strjoin(["euler_", step(i,2), ".jpg"], "")));
%     %close 1 
    figure()
end
k = j;

% Chose the best seperation using RK2 
for j = 1:length(indrk)
    i = indrk(j);
    violation = cell2mat(forplot(i));
    violationrk = cell2mat(forplotrk(i));
    
    plot(cell2mat(lam(i)),log10(violation), pp(j))
    hold on
    plot(cell2mat(lamrk(i)),log10(violationrk), pp(j))
    hold off
    title("Best Total Variation for RK3",'fontsize',16);
    xlabel(strjoin(["Time Step : ", step(i,2)]), 'fontsize',16);
    legend("euler", "rk3")
    %legend(step(indrk,2))
    %hold off
%     savefig(j, char(strjoin(["RK2_", step(i,2), ".fig"], "")));
%     saveas(j,char(strjoin(["RK2_", step(i,2), ".jpg"], "")));
%     close 1 
    figure()
end
 close(k+j+1)
toc()


%% playing around with finding the best time step using euler
n = 20;
x = linspace(0,1,n)';
dx = x(2)-x(1);
nu = 1e-4;

L = nu*D2;  % Linear
N = @(u) (-D1)*u; % Nonlinear or not stiff part

% [violation, lam] = LNTotVar("euler", n,linspace(.01, 2), (-D1+L), @(u) 0*u) % Best case Scenerio
% plot(lam,log10(violation), 'o')
% hold on
[violation, lam] =LNTotVar("euler",n, linspace(.01, 2), 0, @(u) (-D1+L)*u) % Worst case Scenerio
plot(lam,log10(violation), 's')
% LNTotVar("euler", n, linspace(.01, 2), L, @(u) (-D1)*u) % choosing a N value
% plot(lam,log10(violation), '<')
LNTotVar("euler", n, linspace(.01, 2), F, @(u) (-D1+G+H+D+E)*u) % choosing a N value
plot(lam,log10(violation), '<')

%% playing around with finding the best time step using RK2l

[violation, lam] = LNTotVar("rk3", n, linspace(.01, 2), (-D1+L), @(u) 0*u); % Best case Scenerio
plot(lam,log10(violation), 'o')
hold on
[violation, lam] = LNTotVar("rk3", n, linspace(.01, 2), 0, @(u) (-D1+L)*u); % Worst case Scenerio
plot(lam,log10(violation), 'x')
[violation, lam] = LNTotVar("rk3", n, linspace(.01, 2), L, @(u) (-D1)*u); % choosing a N value
plot(lam,log10(violation), 's')
legend('best', 'worst', '-D1 = NonL')