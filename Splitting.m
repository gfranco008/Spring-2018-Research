%--------------------------------------------------------------------------
%% This code determines the best distribution for solving IMEX methods using a Linear/ Nonlinear aprroach
% In this case we are solving u_t = - u_x + nu*u_xx % Here nu is your epsilon

% Gustavo Franco Reynoso
% March 13, 2018
%--------------------------------------------------------------------------

clear all;
close all;

n = 6;
x = linspace(0,1,n)';
dx = x(2)-x(1);
nu = 1;

En1 = diag(ones(n-1,1),-1); % En1 is the subdiagonal of a Downwind matrix
En1(1,end) = 1;

I = diag(ones(n,1));  % Identity matrix

E1 = diag(ones(n-1,1),1); % E1 is the superdiagonal of a Upwind matrix
E1(end,1) = 1;
% when compared to TV_1SRn2ND
% we divided everything even more, so we can use fractions of a matrix
% and see if there is an optimal splitting.

% These are the weights for every diagonal of the 1st order and 2nd order
% FD stencil
w1 = (1/dx); % A = w1 * En1 
w2 = -(1/dx); % B = w2 * I
w3 = (nu*(1/dx^2)); % C = w3 * En1
w4 = (nu*-(2/dx^2)); % D = w4 * I
w5 = (nu*(1/dx^2)); % E = w5 * E1

% We need to save some info for later.
forplotr = [];
forplot1 = [];
steprk = [];
stepr = {};
line = linspace(0,1,n);

tic()

% Looping through all possible combinations of xi
parfor i = 1:n; % has to be integer
    steprk = [];
    stepr{i,:} = {};
    forplotrk{i,:} = {};
    forplot{i,:} = {};
    lines = linspace(0,1,n);
    for j = lines;
        for k = lines;
            for l = lines;
                for m = lines;
                    xi1 = line(i);
                    xi2 = j;
                    xi3 = k;
                    xi4 = l;
                    xi5 = m;
                    % In this case we would put xi(n) in L and [1-xi(n)] in N
                    L = (xi1*w1)*En1 + (xi2*w2)*I + (xi3*w3)*En1 + (xi4*w4)*I + (xi5*w5)*E1;
                    N = (1-xi1)*w1*En1 + (1-xi2)*w2*I + (1-xi3)*w3*En1 + (1-xi4)*w4*I + (1-xi5)*w5*E1;
                    
                    % when the xi(n) equal 1 or 0, it will be the same as TV_1STn2ND
                    % in terms of ABCDE. So they should give me the same value for these exact
                    % points.
                    
                    % The goal is to make it like a sliding bar where you pick 0-1 for every
                    % xi(n)
                    [violationrk, Lambdark] = LNTotVar("rk3", n, linspace(.01, 2), L, @(u) N*u);
                    [violation, Lambda] = LNTotVar("euler", n, linspace(.01, 2), L, @(u) N*u);
                    
                    % This determines the max number of points before there
                    % is a jump in TV
                    pnt = min(find(diff(violationrk)> 1e-4))-1; 
                   % You could also do this for the min amount of non-zeros
                   % That way you could see plots that have jumps in them
                   % early on. 
                    
                   % pnt2 = sum(diff(violationrk)> 1e-4);
                    
                    forplotr = [forplotrk{i,:}, {violationrk}]; % Saves the violation for later use
                    forplot1 = [forplot{i,:}, {violation}];
                    
                    
                   % If pnt is [] or TV is not maintained then it equals zero
                    if isempty(pnt)==1
                        pnt = 0;
                    end
                    % Print out the inputs
                    steprk = [steprk; pnt, string([xi1, xi2, xi3, xi4, xi5])];
                end
            end
        end
    end
    lamrk{i} = Lambdark; % Save Lambda 
    stepr{i,:} = {steprk};
    forplotrk{i,:} = {violationrk}; % Saves the violation for later use
    forplot{i,:} = {violation};
end
toc()

% unpack forplot
forplot = forplot{:};
forplotrk = forplotrk{:};

% Convert cell back to matrix
for a = 1:length(stepr);
    X = (cat(3,stepr{a}));
    steprk = [steprk; X{:}];
end

% Find the indexes for all violations with the max steps
indrk = find(str2double(steprk(:,1)) >= max(str2double(steprk(:,1))));

for j = 1:length(indrk)
    i = indrk(j);
    
    violationrk = forplotrk{i};
    violation = forplot{i};
    
    plot(lamrk,log10(violationrk), 'O')
    hold on
    plot(lamrk,log10(violation))
    title("Best Total Variation for RK3",'fontsize',16);
    % Plots Linear and Non-linear values
    xlabel({"Time Step "; strjoin([" L =", steprk(i,2:6)]);  strjoin([" N =", string(1-str2double(steprk(i,2:6)))])}, 'fontsize',16);
    legend("RK3", "RK2", "Euler")
    pause
    
    %Every 50 graphs, use a new figure.
    if mod(j,50) == 0; 
        figure()
    end
    
end