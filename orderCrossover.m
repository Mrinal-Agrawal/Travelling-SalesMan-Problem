function [c1,c2] = orderCrossover(p1,p2)
n = length(p1);
cp = sort(randperm(n,2));
c1 = zeros(1,n);
c2 = zeros(1,n);
c1(cp(1):cp(2)) = p1(cp(1):cp(2));
c2(cp(1):cp(2)) = p2(cp(1):cp(2));
fill1 = setdiff(p2, c1, 'stable');
fill2 = setdiff(p1, c2, 'stable');
idx = [1:cp(1)-1, cp(2)+1:n];
c1(idx) = fill1; c2(idx) = fill2;
end
