
G = [[  0;   0; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 1;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  1;   1; 0;   0;   0;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   1;   0;   1; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   1; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 0;   0;   1;   0; 0; 0; 0; 0; 0],...
     [  0;   0; 0;   0;   1;   0; 0; 0; 0; 0; 0]];
 
alpha = 0.85;
 
   %Compute the number of pages and store it as the constant R
   R=size(G,2);

   %Create a Rx1 column vector of ones, named e
   e=ones(R,1);
   
   %Create a 1xR row vector, named Deg, to store the degrees of pages
   Deg = sum(G);
   
   %Create a row vector, named D, of the length of the number of dead pages
   %   in the given web to store the indices of the dead pages
   d = find(Deg==0);
   
   %Create a cell, named Links, of R arrays such that each array i stores 
   %    some j's representing the indices of the pages in the web that have
   %    an outlink to the page indexed at i, and their corresponding 
   %    1/deg(j) values. 
   %Each j comes before its 1/deg(j) value, which given an array i in 
   %    Links, j's and 1/deg(j)'s are layed out in a staggered fashipn.
   %Retrieval of the information can be done by picking odd numbered 
   %    elements for indices j's, and even numbered elemnt for 1/deg(j)'s
   %    e.g Links{1}(1:2:end) gives all odd numbered elements in Links{1},
   %        that are the indices j's of those pages that have an outlink to 
   %        page 1; Links{1}(2:2:end) gives all even number elements in
   %        Links{1}, that are the corresponding 1/deg(j) values
   %(initialized to contain R empty matrices)
   Links = cell([1,R]);
   for i = 1:R
      for j = 1:R
	     %If j has an outlink to i, add j's index and 1/deg(j) to C{i}
         if G(i,j)==1
		   Links{i} = [Links{i},j, 1/Deg(j)];
         end
      end
   end
     
   %Initialize the PageRank probablity vector, named p, to e./R
   p = e./R;
   
   %Initialize the iterator to 1
   iter = 1;
   
   while true
       %p_prev preserves p before next iteration 
       p_prev = p;
       
       %pagerank iteration
       p=alpha.*(cellfun(@(x)sum(p(x(1:2:end)).*(x(2:2:end).')),Links).')+...
            alpha.*e.*(sum(p(d))/R)+ ...
            (1-alpha).*e.*(1/R);
        
       %Break out the loop if the PageRank probablity 
       %    vector converges under the tolerance
       if max(abs(p-p_prev)) < 10E-7
               break
       end
       iter = iter + 1;
   end
   
%outputs p and iter for comparison
p
iter
   
%-------------------------------------------------------------------------

%Matrix multiplication approach for comparison purposes
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
 
% create column vector D
D = zeros(11,1);
D(1)=1;

% create matrix Q
Q = P + e*(D.').*(1/11);

% create Google Matrix M (alpha is a variable)
syms alpha
M = (alpha.*Q) + (1-alpha).*(e*(e.')).*(1/11);

% use alpha=0.75
alpha = 0.85;
M=subs(M);

% initialize the ranking vector x with e./(11)
x = e./(11);

%the iter from the O(R) approach is used here for calculation comparison
for i=1:iter
    x = M*x;
end

%Output for comparison
double(x)

   
   
   
   
