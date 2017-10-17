load('math_uwaterloo.mat')

alpha = 0.85;

[p,iter]= MyPageRank(G, alpha);

[y, I] = sort(p,'descend');

disp('Rank obtained with alpha = 0.85');

for n = 1:min(length(I), 20)
    disp([num2str(n) ': ' U{I(n)}]);
end
