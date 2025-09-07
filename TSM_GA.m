% clc; clear; close all;
% 
% %% ----------- Load TSP Dataset -----------
% filename = 'att48.tsp';      % dataset file
% [coords, D] = readTSP(filename);
% nCities = size(coords,1);
% 
% %% ----------- Load Optimal Tour (if available) -----------
% optTourFile = 'att48.opt.tour';   % change to your .tour file name
% if isfile(optTourFile)
%     optTour = readTour(optTourFile);
%     optCost = evaluateTour(optTour, D);
% else
%     optTour = [];
%     optCost = NaN;
% end
% 
% %% ----------- GA Parameters -----------
% popSize = 200;         % number of chromosomes
% maxGen  = 2000;        % number of generations
% pc      = 0.8;        % crossover probability
% pm      = 0.2;        % mutation probability
% 
% %% ----------- Initialize Population -----------
% population = zeros(popSize, nCities);
% for i = 1:popSize
%     population(i,:) = randperm(nCities);
% end
% 
% bestCost = inf;
% bestTour = [];
% 
% %% ----------- Main GA Loop -----------
% for gen = 1:maxGen
%     
%     % ---- Evaluate fitness ----
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
%     % ---- Selection + Crossover + Mutation ----
%     newPop = zeros(popSize, nCities);
%     for i = 1:2:popSize
%         p1 = tournamentSelection(population, fitness);
%         p2 = tournamentSelection(population, fitness);
%         
%         if rand < pc
%             [c1, c2] = orderCrossover(p1, p2);
%         else
%             c1 = p1; c2 = p2;
%         end
%         
%         if rand < pm, c1 = swapMutation(c1); end
%         if rand < pm, c2 = swapMutation(c2); end
%         
%         newPop(i,:)   = c1;
%         newPop(i+1,:) = c2;
%     end
%     population = newPop;
%         
%     % ---- Display progress ----
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
% 
% %% ----------- Functions -----------
% 
% function [coords, distMatrix] = readTSP(filename)
%     fid = fopen(filename, 'r');
%     if fid == -1, error('Cannot open file'); end
%     
%     nCities = 0; coords = [];
%     while true
%         line = fgetl(fid);
%         if contains(line,'DIMENSION')
%             tokens = regexp(line,'\d+','match');
%             nCities = str2double(tokens{1});
%         elseif contains(line,'NODE_COORD_SECTION')
%             coords = zeros(nCities,2);
%             for i=1:nCities
%                 data = fscanf(fid,'%d %f %f',3);
%                 coords(i,:) = data(2:3);
%             end
%             break;
%         end
%     end
%     fclose(fid);
%     
%     distMatrix = zeros(nCities);
%     for i=1:nCities
%         for j=1:nCities
%             if i ~= j
%                 dx = coords(i,1)-coords(j,1);
%                 dy = coords(i,2)-coords(j,2);
% %             if strcmp(edgeType,'EUC_2D')
% %                 distMatrix(i,j) = round(sqrt(dx^2 + dy^2));
% %             elseif strcmp(edgeType,'ATT')
%                 rij = sqrt((dx^2 + dy^2)/10.0);
%                 tij = round(rij);
%                 if tij < rij
%                     distMatrix(i,j) = tij + 1;
%                 else
%                     distMatrix(i,j) = tij;
%                 end
% %             end
%                % distMatrix(i,j) = round(sqrt(dx^2 + dy^2));
%             end
%         end
%     end
% end
% 
% 
% 
% 
% function cost = evaluateTour(tour, D)
%     n = length(tour);
%     cost = 0;
%     for i=1:n-1
%         cost = cost + D(tour(i), tour(i+1));
%     end
%     cost = cost + D(tour(end), tour(1));
% end
% 
% function p = tournamentSelection(pop, fitness)
%     k = 3;
%     idx = randperm(size(pop,1),k);
%     [~,bestIdx] = min(fitness(idx));
%     p = pop(idx(bestIdx),:);
% end
% 
% function [c1,c2] = orderCrossover(p1,p2)
%     n = length(p1);
%     cp = sort(randperm(n,2));
%     c1 = zeros(1,n); c2 = zeros(1,n);
%     c1(cp(1):cp(2)) = p1(cp(1):cp(2));
%     c2(cp(1):cp(2)) = p2(cp(1):cp(2));
%     fill1 = setdiff(p2, c1, 'stable');
%     fill2 = setdiff(p1, c2, 'stable');
%     idx = [1:cp(1)-1, cp(2)+1:n];
%     c1(idx) = fill1; c2(idx) = fill2;
% end
% 
% function c = swapMutation(p)
%     n = length(p);
%     idx = randperm(n,2);
%     c = p;
%     c(idx) = c(fliplr(idx));
% end
% 
% function tour = readTour(filename)
%     fid = fopen(filename,'r');
%     if fid == -1, error('Cannot open file'); end
%     
%     tour = [];
%     while true
%         line = fgetl(fid);
%         if contains(line,'TOUR_SECTION')
%             break;
%         end
%     end
%     
%     while true
%         line = fgetl(fid);
%         val = str2double(line);
%         if val == -1 || isnan(val), break; end
%         tour(end+1) = val; %#ok<AGROW>
%     end
%     fclose(fid);
% end


clc; clear; close all;

%% ----------- Load TSP Dataset -----------
filename = 'a280.tsp';      
[coords, D] = readTSP(filename);
nCities = size(coords,1);

%% ----------- Load Optimal Tour (if available) -----------
optTourFile = 'a280.opt.tour';   
if isfile(optTourFile)
    optTour = readTour(optTourFile);
    optCost = evaluateTour(optTour, D);
else
    optTour = [];
    optCost = NaN;
end

%% ----------- GA Parameters -----------
popSize = 300;         % Increased population
maxGen  = 500;        % More generations
pc      = 0.8;         
pm      = 0.3;         % Slightly higher mutation
tournamentK = 2;       % Stronger selection pressure
elitismN = 1;          % Number of elites to keep

%% ----------- Initialize Population -----------
population = zeros(popSize, nCities);
for i = 1:popSize
    population(i,:) = randperm(nCities);
end

bestCost = inf;
bestTour = [];

%% ----------- Main GA Loop -----------
for gen = 1:maxGen
    
    % ---- Evaluate fitness ----
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
    
    % ---- Selection + Crossover + Mutation + Local Search ----
    newPop = zeros(popSize, nCities);
    
    % ---- Elitism ----
    elites = population(idx,:);
    newPop(1:elitismN,:) = elites;
    
    for i = elitismN+1:2:popSize
        p1 = tournamentSelection(population, fitness, tournamentK);
        p2 = tournamentSelection(population, fitness, tournamentK);
        
        if rand < pc
            [c1, c2] = orderCrossover(p1, p2);
        else
            c1 = p1; c2 = p2;
        end
        
        if rand < pm, c1 = inversionMutation(c1); end
        if rand < pm, c2 = inversionMutation(c2); end
        
        % ---- 2-opt local search ----
        c1 = twoOpt(c1, D);
        c2 = twoOpt(c2, D);
        
        newPop(i,:) = c1;
        if i+1 <= popSize
            newPop(i+1,:) = c2;
        end
    end
    population = newPop;
    
    % ---- Display progress ----
    if ~isnan(optCost)
        fprintf('Gen %d | Best = %d | Optimum = %d | Gap = %.2f%%\n', ...
            gen, bestCost, optCost, 100*(bestCost-optCost)/optCost);
    else
        fprintf('Gen %d | Best = %d\n', gen, bestCost);
    end
end

%% ----------- Plot Best Tour -----------
figure;
plot(coords(bestTour,1), coords(bestTour,2), 'b-o','LineWidth',1.5);
hold on;
plot([coords(bestTour(end),1) coords(bestTour(1),1)], ...
     [coords(bestTour(end),2) coords(bestTour(1),2)], 'b-o','LineWidth',1.5);

if ~isempty(optTour)
    plot(coords(optTour,1), coords(optTour,2), 'r--','LineWidth',1.5);
    plot([coords(optTour(end),1) coords(optTour(1),1)], ...
         [coords(optTour(end),2) coords(optTour(1),2)], 'r--','LineWidth',1.5);
    legend('GA Best Tour',['Optimal Tour (Cost = ' num2str(optCost) ')']);
else
    legend('GA Best Tour');
end

title(['Best GA Tour (Cost = ', num2str(bestCost), ')']);
xlabel('X'); ylabel('Y'); grid on;

%% ----------- Functions -----------

function [coords, distMatrix] = readTSP(filename)
    fid = fopen(filename, 'r');
    if fid == -1, error('Cannot open file'); end
    nCities = 0; coords = [];
    while true
        line = fgetl(fid);
        if contains(line,'DIMENSION')
            tokens = regexp(line,'\d+','match');
            nCities = str2double(tokens{1});
        elseif contains(line,'NODE_COORD_SECTION')
            coords = zeros(nCities,2);
            for i=1:nCities
                data = fscanf(fid,'%d %f %f',3);
                coords(i,:) = data(2:3);
            end
            break;
        end
    end
    fclose(fid);
    
    distMatrix = zeros(nCities);
    for i=1:nCities
        for j=1:nCities
            if i ~= j
                dx = coords(i,1)-coords(j,1);
                dy = coords(i,2)-coords(j,2);
%                 rij = sqrt((dx^2 + dy^2)/10.0);
%                 tij = round(rij);
%                 if tij < rij
%                     distMatrix(i,j) = tij + 1;
%                 else
%                     distMatrix(i,j) = tij;
%                 end
                 distMatrix(i,j) = round(sqrt(dx^2 + dy^2));
            end
        end
    end
end

function cost = evaluateTour(tour, D)
    n = length(tour);
    cost = 0;
    for i=1:n-1
        cost = cost + D(tour(i), tour(i+1));
    end
    cost = cost + D(tour(end), tour(1));
end

function p = tournamentSelection(pop, fitness, k)
    idx = randperm(size(pop,1),k);
    [~,bestIdx] = min(fitness(idx));
    p = pop(idx(bestIdx),:);
end

function [c1,c2] = orderCrossover(p1,p2)
    n = length(p1);
    cp = sort(randperm(n,2));
    c1 = zeros(1,n); c2 = zeros(1,n);
    c1(cp(1):cp(2)) = p1(cp(1):cp(2));
    c2(cp(1):cp(2)) = p2(cp(1):cp(2));
    fill1 = setdiff(p2, c1, 'stable');
    fill2 = setdiff(p1, c2, 'stable');
    idx = [1:cp(1)-1, cp(2)+1:n];
    c1(idx) = fill1; c2(idx) = fill2;
end

function c = inversionMutation(p)
    n = length(p);
    idx = sort(randperm(n,2));
    c = p;
    c(idx(1):idx(2)) = fliplr(c(idx(1):idx(2)));
end

function tour = twoOpt(tour, D)
    n = length(tour);
    improved = true;
    while improved
        improved = false;
        for i = 1:n-2
            for j = i+2:n
                if j == n && i == 1, continue; end
                a = tour(i); b = tour(i+1);
                c = tour(j); d = tour(mod(j,n)+1);
                delta = (D(a,c)+D(b,d)) - (D(a,b)+D(c,d));
                if delta < 0
                    tour(i+1:j) = fliplr(tour(i+1:j));
                    improved = true;
                end
            end
        end
    end
end

function tour = readTour(filename)
    fid = fopen(filename,'r');
    if fid == -1, error('Cannot open file'); end
    tour = [];
    while true
        line = fgetl(fid);
        if contains(line,'TOUR_SECTION'), break; end
    end
    while true
        line = fgetl(fid);
        val = str2double(line);
        if val == -1 || isnan(val), break; end
        tour(end+1) = val; %#ok<AGROW>
    end
    fclose(fid);
end

