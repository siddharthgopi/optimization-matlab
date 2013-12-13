clear;

gridsize = input('Dr.Alex, please specify a size for the k-grid; ')
eps = input('Please specify your precision parameter; ')

% parameters
beta = 0.95 ;
gamma = 3.2938;
delta = 0.1;
theta = 0.33;
A= 1 ;

%calculate steady state
kstar=((1/beta-(1-delta))/(A*theta))^(1/(theta-1))
cstar = A*kstar^theta-delta*kstar

%intialize k-grid
N=gridsize;

klo=kstar*0.9;
khi=kstar*1.1; 

step =(khi-klo)/N;
k = klo:step:khi;
n=length(k);

%creating a matrix for V0 with proper dimensions
%first guess assumes all k is immediately consumed
ktheta =k.^theta;
colones = ones(n,1);  %create a column vector of ones
s = colones*ktheta;   %create a matrix of kthetas 
s1= colones*k;        %create another matrix of k's

ytot = s'+(1-delta)*s1';  % the resultant matrix here represents -- 
                          % k^theta + (1-d)k - kt+1. But
                          % kt+1 is 0 since our first guess.

v =(ytot.^(1-gamma)-1)/(1-gamma); %applying the u function to this matrix

%create a part of the matrix that does not depend on Bw(kt+1)
rowones = colones';
I = k'*rowones;
J = colones*k;
C = (J.^theta) + (1-delta)*J - I;  %nxn matrix

C2= (k'.^theta) + (1-delta)*k';    %this particular vector was created to
                                   %help calculate vector representing
                                   %stable arm for later use.
                                   
U = (C.^(1-gamma)-1)/(1-gamma);  % util levels for each consumption level


r = U+beta*v;    % adding the time variant part back
v1 = max(r);     % Row vector of only the maximum of each column of r

%The main loop
change = eps;

while change >= eps             % while distance is above a threshold
    v1old = v1;                 % assigning v1old to calculate euc.distance later
    w= ones(n,1)*v1;            % since v1 is a row vector, we make it a matrix
    w1 = U+beta*w';             % adding value of carrying forward kt+1_i level of capital
    v1 = max(w1);               % getting back the v1 vector for next iteration
    change = norm(v1-v1old,2);  % Finds euclidean distance of vectors
    
    change;                     %i turn this on when im playing around with 
                                % parameter values to see the convergence
    
end

[val,ind] = max(w1); %save the values and corresponding index

optk = k(ind);       %this represents the optimal capital accumulation policy

optc = C2-optk';     %optimal consumpion values, applied the formula from the 
                     %langrangian here, transposing optk to get correct
                     %dimensions
                     

deltak0 = A*k.^theta-delta*k; %represents curve ct = f(kt) - dkt, i.e. change in k is 0

deltac0=deltak0+k-kstar;      %this represents points where,
                              %change in consumption is 0

                              

%**** Plotting g(k), Stable Arm *****%

wish = input('Press ''y'' to view Optimal Policy for Capital Accumulation, g(k), and the Stable Arm ','s');

if wish == 'y'
    
figure(1);    
plot(k,optk',k,k,'linewidth',1); hold on;
plot(kstar, kstar, 'r+');                        % put a red plus on the solution point
xlabel('Current capital')
ylabel('Optimal future capital')
text(kstar+.1,kstar-.1,['k*=' num2str(kstar)]);  % text with solution offset by 0.1 to the right and up
title('Optimal Policy function, g(k)')                                          

figure(2);
plot(k,optc,k,deltak0,k,deltac0)
xlabel('k_t')
ylabel('c_t')
legend('Stable Arm','Delta k_t = 0','Delta c_t = 0') 
title('Phase Diagram')

figure(3);
plot(k,optc)
xlabel('k_t')
ylabel('c_t')
legend('Stable Arm')
title('Stable Arm Singled out')


figure(4);
plot(k,v1)
xlabel('Choice of kt+1')
ylabel('Value function')
title('Value function against k')


else disp('End')
 
end




