clear;

iter = input('Dr Alex, please specify number of iterations ')
gridsize = input('Please specify the size of the k-grid ')

% parameters
beta = 0.96 ;
gamma = 2 ;
delta = 0.05;
theta = 0.34;
A= 1 ;

%calculate steady state
kstar=((1/beta-(1-delta))/(A*theta))^(1/(theta-1));

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

ytot = s'+(1-delta)*s1';  %this is matrix of -- ktheta + (1-d)k -- minus kt+1 but it is 0

v =(ytot.^(1-gamma)-1)/(1-gamma); %applying the u function to this matrix

%create a part of the matrix that does not depend on Bw(kt+1)
rowones = colones';
I = k'*rowones;
J = colones*k;
C = (J.^theta) + (1-delta)*J - I;

U = (C.^(1-gamma)-1)/(1-gamma);  % u levels for each c level

%rule for choosing next iter
r = U+beta*v;
v1 = max(r);


%actual iterations
for j=1:iter
    w=ones(n,1)*v1;
    w1 = U+beta*w';
    v1 = max(w1);
   
end

[val,ind] = max(w1); %save the values and corresponding index
optk = k(ind);   

wish = input('Press ''y'' to view plot ','s')
%plot k against optimal k'
if wish == 'y'
    
plot(k,optk',k,k,':','linewidth',2)
xlabel('Current capital')
ylabel('Optimal future capital')


else
end




