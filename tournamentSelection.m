function p = tournamentSelection(pop, fitness, k)
idx = randperm(size(pop,1),k);
[~,bestIdx] = min(fitness(idx));
p = pop(idx(bestIdx),:);
end
