load('math_uwaterloo.mat')

Alpha = [0.15,0.6,0.75,0.95];
times = size(Alpha,2);

for i=1: times
    fprintf('Rank obtained with alpha = %g\n', Alpha(i));
    [p,iter]= MyPageRank(G, Alpha(i));
    [y, I] = sort(p,'descend');
    for n = 1:min(length(I), 20)
        disp([num2str(n) ': ' U{I(n)}]);
    end
    fprintf('which takes %d iterations\n',iter);
end
