clc; clear; close all;

%% ----------- Load TSP Dataset -----------
filename = 'att48.tsp';      
[coords, D] = readTSP(filename);
nCities = size(coords,1);

%% ----------- Load Optimal Tour (if available) -----------
optTourFile = 'att48.opt.tour';   
if isfile(optTourFile)
    optTour = readTour(optTourFile);
    optCost = evaluateTour(optTour, D);
else
    optTour = [];
    optCost = NaN;
end

% %% ----------- GA Parameters -----------
popSize = 400;
maxGen  = 3000;
pc      = 0.8;
pm      = 0.3;
tournamentK = 2;
elitismN = 1;
% 
%% ----------- Initialize Population -----------
population = zeros(popSize, nCities);
for i = 1:popSize
    population(i,:) = randperm(nCities);
end

bestCost = inf;
bestTour = [];

% %% ----------- Main GA Loop -----------
% for gen = 1:maxGen
%     
%     % Evaluate fitness
%     fitness = zeros(popSize,1);
%     for i = 1:popSize
%         fitness(i) = evaluateTour(population(i,:), D);
%     end
%     
%     % Track best solution
%     [minCost, idx] = min(fitness);
%     if minCost < bestCost
%         bestCost = minCost;
%         bestTour = population(idx,:);
%     end
%     
%     % Selection, Crossover, Mutation, Local Search
%     newPop = zeros(popSize, nCities);
%     
%     % Elitism
%     elites = population(idx,:);
%     newPop(1:elitismN,:) = elites;
%     
%     for i = elitismN+1:2:popSize
%         p1 = tournamentSelection(population, fitness, tournamentK);
%         p2 = tournamentSelection(population, fitness, tournamentK);
%         
%         if rand < pc
%             [c1, c2] = orderCrossover(p1, p2);
%         else
%             c1 = p1; c2 = p2;
%         end
%         
%         if rand < pm, c1 = inversionMutation(c1); end
%         if rand < pm, c2 = inversionMutation(c2); end
%         
%         % 2-opt local search
%         c1 = twoOpt(c1, D);
%         c2 = twoOpt(c2, D);
%         
%         newPop(i,:) = c1;
%         if i+1 <= popSize
%             newPop(i+1,:) = c2;
%         end
%     end
%     population = newPop;
%     
%     % Display progress
%     if ~isnan(optCost)
%         fprintf('Gen %d | Best = %d | Optimum = %d | Gap = %.2f%%\n', ...
%             gen, bestCost, optCost, 100*(bestCost-optCost)/optCost);
%     else
%         fprintf('Gen %d | Best = %d\n', gen, bestCost);
%     end
% end
% 
% %% ----------- Plot Best Tour -----------
% figure;
% plot(coords(bestTour,1), coords(bestTour,2), 'b-o','LineWidth',1.5);
% hold on;
% plot([coords(bestTour(end),1) coords(bestTour(1),1)], ...
%      [coords(bestTour(end),2) coords(bestTour(1),2)], 'b-o','LineWidth',1.5);
% 
% if ~isempty(optTour)
%     plot(coords(optTour,1), coords(optTour,2), 'r--','LineWidth',1.5);
%     plot([coords(optTour(end),1) coords(optTour(1),1)], ...
%          [coords(optTour(end),2) coords(optTour(1),2)], 'r--','LineWidth',1.5);
%     legend('GA Best Tour',['Optimal Tour (Cost = ' num2str(optCost) ')']);
% else
%     legend('GA Best Tour');
% end
% 
% title(['Best GA Tour (Cost = ', num2str(bestCost), ')']);
% xlabel('X'); ylabel('Y'); grid on;




%% ----------- Main GA Loop (μ + λ Selection) -----------
for gen = 1:maxGen
    
    % ---- Evaluate fitness of current population ----
    fitness = zeros(popSize,1);
    for i = 1:popSize
        fitness(i) = evaluateTour(population(i,:), D);
    end
    
    % Track best solution
    [minCost, idx] = min(fitness);
    if minCost < bestCost
        bestCost = minCost;
        bestTour = population(idx,:);
    end


    
    % ---- Generate offspring ----
    newPop = zeros(popSize, nCities);
    for i = 1:2:popSize
        p1 = tournamentSelection(population, fitness, 3);
        p2 = tournamentSelection(population, fitness, 3);
        
        if rand < pc
            [c1, c2] = orderCrossover(p1, p2);
        else
            c1 = p1; c2 = p2;
        end
        
        if rand < pm, c1 = inversionMutation(c1); end
        if rand < pm, c2 = inversionMutation(c2); end
        
%         optional local search improvement
        c1 = twoOpt(c1, D);
        c2 = twoOpt(c2, D);
        
        newPop(i,:) = c1;
        if i+1 <= popSize
            newPop(i+1,:) = c2;
        end
    end
    
    % ---- Combine old + new populations (μ + λ) ----
    combinedPop = [population; newPop];              % 2n × nCities
    fitnessCombined = zeros(2*popSize,1);
    for i = 1:2*popSize
        fitnessCombined(i) = evaluateTour(combinedPop(i,:), D);
    end
    
    % ---- Select the best n ----
    [~, idxSorted] = sort(fitnessCombined);
    population = combinedPop(idxSorted(1:popSize), :);
    
    % ---- Display progress ----
    if ~isnan(optCost)
        fprintf('Gen %d | Best = %d | Optimum = %d | Gap = %.2f%%\n', ...
            gen, bestCost, optCost, 100*(bestCost-optCost)/optCost);
    else
        fprintf('Gen %d | Best = %d\n', gen, bestCost);
    end
end

