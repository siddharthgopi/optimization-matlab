clear;


gridsize = input('Please specify a size for the k-grid; ')
eps = input('Please specify your precision parameter; ')


% parameters
beta = 0.96 ;
gamma = 2 ;
delta = 0.05;
theta = 0.34;
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
C = (J.^theta) + (1-delta)*J - I;
C2= (k'.^theta) + (1-delta)*k';

U = (C.^(1-gamma)-1)/(1-gamma);  % u levels for each c level

%rule for choosing next iter
r = U+beta*v;
v1 = max(r);

%actual iterations
change = eps;

while change >= eps
    v1old = v1;
    w= ones(n,1)*v1;
    w1 = U+beta*w';
    v1 = max(w1);
    change = norm(v1-v1old,2);
    
    change;
    
   
end

[val,ind] = max(w1); %save the values and corresponding index
optk = k(ind);       %this represents the optimal capital accumulation policy

optc = C2-optk';     %optimal consumpion values, applied the formula from the 
                     %langrangian here, transposing optk to get correct
                     %dimensions
                     

deltak0 = A*k.^theta-delta*k; %represents curve ct = f(kt) - dkt
deltac0=deltak0+k-kstar;      %similarly this represents points where ...
                              %change in consumption is 0

                              
%plot k against optimal k'

wish = input('Press ''y'' to view Optimal Policy for Capital Accumulation, g(k), and the Stable Arm ','s');
if wish == 'y'
    
figure(1);    
plot(k,optk',k,k,'linewidth',1); hold on;
xlabel('Current capital')
ylabel('Optimal future capital')

figure(2);
plot(k,optc,k,deltak0,k,deltac0)


else disp('End')
    
end




