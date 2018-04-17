  clear all;
  close all;

n = 10;
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
w4 = (nu*(-2/dx^2)); % D = w4 * I
w5 = (nu*(1/dx^2)); % E = w5 * E1

% We need to save some info for later.
forplotrk = {};
forplot = {};
steprk = [];
stepr = {};
line = linspace(0,1,n);

xi = [0, 1];
z = 1;

for i = xi
    for j = xi
        for k = xi
            for l = xi
                for m = xi
                    
                    xi1 = i;
                    xi2 = j;
                    xi3 = k;
                    xi4 = l;
                    xi5 = m;
                    
                    L = (xi1*w1)*En1 + (xi2*w2)*I + (xi3*w3)*En1 + (xi4*w4)*I + (xi5*w5)*E1;
                    N = (1-xi1)*w1*En1 + (1-xi2)*w2*I + (1-xi3)*w3*En1 + (1-xi4)*w4*I + (1-xi5)*w5*E1;
                    
                    [violationrk, Lambdark] = LNTotVar("rk3", n, linspace(.01, 4), L, @(u) N*u);
                    
                    pnt = min(find(diff(violationrk)> 1e-4))-1;
                    
                    forplotrk{z} = violationrk; % Saves the violation for later use
                    
                    % If pnt is [] or TV is not maintained then it equals zero
                    if isempty(pnt)==1
                        pnt = 0;
                    end
                    % Print out the inputs
                    steprk = [steprk; pnt, string([xi1, xi2, xi3, xi4, xi5]),z];
                    
%                     figure()
%                     plot(Lambdark,log10(violationrk), 'O')
%                     title("Total Variation for RK3",'fontsize',16);
%                     % Plots Linear and Non-linear values
%                     xlabel({"Time Step "; strjoin([" L =", steprk(z,2:6)]);  strjoin([" N =", string(1-str2double(steprk(z,2:6)))])}, 'fontsize',16);
%                     legend("RK3")
                    
                    %figure()
                    
                    z = z + 1;
                end
            end
        end
    end
    
    
    
end



ind = [1, 32, 9, 3, 27, 16, 8, 25, 11]  
    
legends = ["L = 0"
"L = D1 + D2"
"L = diag(D1)"
"L = diag(D2)"
"L = D1 + diag(D2)"
"L = diag(D1) + D2"
"L = D2"
"L = D1"
"L = diag(D1) + diag(D2)"
"L = IMEX"
]
%figure()
z = 1;
mark = ['s' 'o' '+' ':*' 'd' '<' 'h' '-.' 'x' '-p' '-h'];
for i = ind 
    violationrk = cell2mat(forplotrk(i));
    subplot(3,4,z)
    plot(Lambdark,log10(violationrk), mark(z))
                    title(legends(z),'fontsize',16);
                    % Plots Linear and Non-linear values
                    xlabel({"Time Step "; strjoin([" L =", steprk(i,2:6)]);  strjoin([" N =", string(1-str2double(steprk(i,2:6)))])}, 'fontsize',16);
                    %hold on
                    strjoin([" L =", steprk(i,2:6)])  
                    strjoin([" N =", string(1-str2double(steprk(i,2:6)))])
                    
                
    z = z+1;
end
% legend("L = 0", "L = D1 + D2", "L = diag(D1)", "L = diag(D2)", "L = diag(D1) + diag(D2)", "L = diag(D1) + D2", "L = D1 + diag(D2)", "L = D1", "L = D2")
                    