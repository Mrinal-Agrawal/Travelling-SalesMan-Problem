function c = inversionMutation(p)
n = length(p);
idx = sort(randperm(n,2));
c = p;
c(idx(1):idx(2)) = fliplr(c(idx(1):idx(2)));
end
