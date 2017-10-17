% create P matrix
P = [[  0;   0; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 1;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [1/2; 1/2; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0; 1/3; 0; 1/3;   0; 1/3; 0; 0; 0; 0; 0],...
     [  0; 1/2; 0;   0; 1/2;   0; 0; 0; 0; 0; 0],...
     [  0; 1/2; 0;   0; 1/2;   0; 0; 0; 0; 0; 0],...
     [  0; 1/2; 0;   0; 1/2;   0; 0; 0; 0; 0; 0],...
     [  0; 1/2; 0;   0; 1/2;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 0;   0;   1;   0; 0; 0; 0; 0; 0]];
 
% create e column vector
e = ones(11,1);

% create d column vector
d = zeros(11,1);
d(1)=1;

% create Q matrix
Q = P + e*(d.').*(1/11);

% create Google Matrix M (alpha is a variable)
syms alpha
M = (alpha.*Q) + (1-alpha).*(e*(e.')).*(1/11);

% use alpha=0.75
alpha = 0.75;
M=subs(M);

% initialize the probablity vector p with e./(11)
p = e./(11);

% set interation time = 15
iter = 15;


for i = 1:iter
   p = M*p;
end

disp('With alpha=0.75 and 15 iterations, the ranking vector x =');
disp(double(p));
