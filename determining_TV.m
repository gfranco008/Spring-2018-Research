n = 40;
x = linspace(0,1,n)';
dx = x(2)-x(1);
nu = 1e-4;

% D1= diag(ones(n,1))+ diag(-ones(n-1,1),-1);
% D1(1,end)=-1;
% D1=(1/dx)*D1;

% -- 3 points 
D1 = 2*diag(ones(n,1))+ diag(-(3/2)*ones(n-1,1),-1) + diag(-(1/2)*ones(n-1,1),1);
D1(1, end) = -3/2;
D1(end, 1) = -1/2;
D1 = (1/dx)*D1;


% %--second-order derivative
% D2 = -2*diag(ones(n,1))+ diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
% D2(end,1) = 1;
% D2(1,end) = 1;
% D2 = (1/dx)^2*D2;

% - -- - Centered diff with 4 pnts

D2 = -(5/2)*diag(ones(n,1)) + (4/3)*diag(ones(n-1,1),1) + (4/3)*diag(ones(n-1,1),-1) + -(1/12)*diag(ones(n-2,1),-2) + -(1/12)*diag(ones(n-2,1),2);
D2(end,1) = (4/3);
D2(1,end) = (4/3);
D2(1,end-1) = -(1/12);
D2(2,end) = -(1/12);
D2(end-1,1) = -(1/12);
D2(end,2) = -(1/12);
D2 = nu*(1/dx)^2*D2;
% Split into linear and 'nonlinear'(or whatever you are not treating exactly)

A = diag(-(3/2)*ones(n-1,1),-1);
A(1,end) = -(3/2);
A = -(1/dx)*A;

B = 2*diag(ones(n,1));
B = -(1/dx)*B;

C = diag(-(1/2)*ones(n-1,1),1);
C(end,1) = -(1/2);
C = -(1/dx)*C;


D = -(1/12)*diag(ones(n-2,1),-2);
D(1,end-1) = -(1/12);
D(2,end) = -(1/12);
D = nu*(1/dx)^2*D;

E = (4/3)*diag(ones(n-1,1),-1);
E(1,end) = 4/3;
E = nu*(1/dx)^2*E;

F = -(5/2)*diag(ones(n,1));
F = nu*(1/dx)^2*F;


G = (4/3)*diag(ones(n-1,1),1);
G(end, 1) = 4/3;
G = nu*(1/dx)^2*G;


H = -(1/12)*diag(ones(n-2,1),2);
H(end-1,1) = -(1/12);
H(end,2) = -(1/12);
H = nu*(1/dx)^2*H;

fbp = {A, B, C, D, E, F, G ,H};
step = [];
steprk = [];
forplot = {};
forplotrk = {};
m = length(fbp);
k = 0;
lam = {};
lamrk = {};

for i = 1:length(fbp)-1
    comb = nchoosek(1:m,i);
    j = 0;
    while j < length(comb)
        j = j+1;
        k = k+1;
        comb(j,:);
        fbp = {A, B, C, D, E, F, G ,H};
        ofbp = fbp(comb(j,:));
        linmat = sum(cat(3,ofbp{:}),3);
        fbp(comb(j,:)) = [];
        nonmat = sum(cat(3,fbp{:}),3);
        [violation, Lambda] = LNTotVar("euler", n, linspace(.01, 2), linmat, @(u) nonmat*u);
        pnt = min(find(diff(violation)> 1e-4))-1;
        forplot{k} = violation;
        lam{k} = Lambda;
        if isempty(pnt)==1
            pnt = 0;
        end
        step = [step; pnt, string(char(comb(j,:) + 64))];
        % hold on
        [violationrk, Lambdark] = LNTotVar("rk2", n, linspace(.01, 2), linmat, @(u) nonmat*u);
        pnt = min(find(diff(violationrk)> 1e-4))-1;
        forplotrk{k} =  violationrk;
        lamrk{k} = Lambdark;
        if isempty(pnt)==1
            pnt = 0;
        end
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
    title("Best Total Variation for RK2",'fontsize',16);
    xlabel(strjoin(["Time Step : ", step(i,2)]), 'fontsize',16);
    legend("euler", "rk2")
    %legend(step(indrk,2))
    %hold off
%     savefig(j, char(strjoin(["RK2_", step(i,2), ".fig"], "")));
%     saveas(j,char(strjoin(["RK2_", step(i,2), ".jpg"], "")));
%     close 1 
    figure()
end
 %close(k+j+1)

%% playing around with finding the best time step using euler
n = 10;
x = linspace(0,1,n)';
dx = x(2)-x(1);
nu = 1e-4;

L = D2;  % Linear
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

[violation, lam] = LNTotVar("rk2", n, linspace(.01, 2), (-D1+L), @(u) 0*u) % Best case Scenerio
plot(lam,log10(violation), 'o')
hold on
[violation, lam] = LNTotVar("rk2", n, linspace(.01, 2), 0, @(u) (-D1+L)*u) % Worst case Scenerio
plot(lam,log10(violation), 'o')
[violation, lam] = LNTotVar("rk2", n, linspace(.01, 2), L, @(u) (-D1)*u) % choosing a N value
plot(lam,log10(violation), 'o')