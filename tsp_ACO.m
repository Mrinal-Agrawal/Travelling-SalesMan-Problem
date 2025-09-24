% % function [bestTour, bestCost] = tsp_ACO(filename, optTourFile, optVal)
% % % Ant Colony Optimization for TSP (TSPLIB EUC_2D)
% % % filename    : TSPLIB .tsp file
% % % optTourFile : (optional) .opt.tour file
% % % optVal      : (optional) known optimum value (if no tour file)
% % %
% % % bestTour    : best tour found
% % % bestCost    : cost of best tour
% % 
% % %% Parameters
% % nAnts  = 30;    % number of ants
% % nIter  = 1000;  % iterations
% % alpha  = 1;     % pheromone influence
% % beta   = 5;     % heuristic influence
% % rho    = 0.5;   % evaporation rate
% % Q      = 100;   % pheromone deposit factor
% % 
% % %% 1. Read coordinates
% % coords = readTSP(filename);
% % n = size(coords,1);
% % 
% % %% 2. Build distance matrix
% % D = zeros(n);
% % for i = 1:n
% %     for j = i+1:n
% %         dx = coords(i,1) - coords(j,1);
% %         dy = coords(i,2) - coords(j,2);
% %         d = sqrt(dx^2 + dy^2);
% %         D(i,j) = floor(d+0.5); 
% %         D(j,i) = D(i,j);
% %     end
% % end
% % 
% % %% 3. Initialize pheromones
% % pheromone = ones(n);
% % 
% % %% 4. ACO main loop
% % bestCost = inf;
% % bestTour = [];
% % 
% % fprintf('Starting ACO on %s with %d ants and %d iterations...\n', filename, nAnts, nIter);
% % for it = 1:nIter
% %     allTours = zeros(nAnts, n+1);
% %     allCosts = zeros(1, nAnts);
% % 
% %     % construct tours
% %     for k = 1:nAnts
% %         allTours(k,:) = constructTour(n,D,pheromone,alpha,beta);
% %         allCosts(k)   = tourLength(allTours(k,:),D);
% %     end
% % 
% %     % find best in this iteration
% %     [minCost, idx] = min(allCosts);
% %     if minCost < bestCost
% %         bestCost = minCost;
% %         bestTour = allTours(idx,:);
% %     end
% % 
% %     % pheromone update
% %     pheromone = (1-rho)*pheromone;
% %     for k = 1:nAnts
% %         Lk = allCosts(k);
% %         tour = allTours(k,:);
% %         for i = 1:n
% %             u = tour(i); v = tour(i+1);
% %             pheromone(u,v) = pheromone(u,v) + Q/Lk;
% %             pheromone(v,u) = pheromone(u,v);
% %         end
% %     end
% % 
% %     % log iteration
% %     fprintf('Iter %4d | Best Cost = %d\n', it, bestCost);
% % end
% % 
% % %% 5. Handle optimal value / tour file
% % if exist('optTourFile','var') && ~isempty(optTourFile) && isfile(optTourFile)
% %     optTour = readOptTour(optTourFile);
% %     optCost = tourLength(optTour,D);
% %     fprintf('\nACO best cost: %d | Known OPT from tour file: %d\n', bestCost, optCost);
% % elseif exist('optVal','var') && ~isempty(optVal)
% %     fprintf('\nACO best cost: %d | Provided OPT value: %d\n', bestCost, optVal);
% % else
% %     fprintf('\nACO best cost: %d (no OPT reference available)\n', bestCost);
% % end
% % 
% % end
% % 
% % %% ---------------- Helper Functions ----------------
% % 
% % function coords = readTSP(filename)
% % fid = fopen(filename);
% % if fid < 0, error('Cannot open file %s',filename); end
% % coords = [];
% % while true
% %     tline = fgetl(fid);
% %     if ~ischar(tline), break; end
% %     if contains(tline,'NODE_COORD_SECTION'), break; end
% % end
% % while true
% %     tline = fgetl(fid);
% %     if ~ischar(tline) || contains(tline,'EOF'), break; end
% %     nums = sscanf(tline,'%d %f %f');
% %     coords(end+1,:) = nums(2:3)'; %#ok<AGROW>
% % end
% % fclose(fid);
% % end
% % 
% % function L = tourLength(tour,D)
% % L = 0;
% % for i = 1:length(tour)-1
% %     L = L + D(tour(i),tour(i+1));
% % end
% % end
% % 
% % function tour = constructTour(n,D,pheromone,alpha,beta)
% % tour = zeros(1,n+1);
% % start = randi(n);
% % tour(1) = start;
% % visited = false(1,n);
% % visited(start) = true;
% % 
% % for i = 2:n
% %     u = tour(i-1);
% %     probs = zeros(1,n);
% %     for v = 1:n
% %         if ~visited(v)
% %             tau = pheromone(u,v)^alpha;
% %             eta = (1/D(u,v))^beta;
% %             probs(v) = tau * eta;
% %         end
% %     end
% %     
% %     if sum(probs) == 0
% %         unvisited = find(~visited);
% %         v = unvisited(randi(numel(unvisited)));
% %     else
% %         probs = probs / sum(probs);
% %         v = rouletteWheel(probs);
% %     end
% % 
% %     tour(i) = v;
% %     visited(v) = true;
% % end
% % 
% % tour(end) = start; % close the loop
% % end
% % 
% % function v = rouletteWheel(probs)
% % probs = probs(:)'; % row
% % if all(probs == 0)
% %     v = randi(length(probs)); % fallback random
% %     return;
% % end
% % r = rand * sum(probs); % scale
% % c = cumsum(probs);
% % v = find(r <= c,1,'first');
% % if isempty(v)
% %     v = length(probs); % last index fallback
% % end
% % end
% % 
% % function tour = readOptTour(filename)
% % fid = fopen(filename);
% % if fid < 0, error('Cannot open file %s',filename); end
% % tour = [];
% % while true
% %     tline = fgetl(fid);
% %     if ~ischar(tline), break; end
% %     if contains(tline,'TOUR_SECTION')
% %         while true
% %             tline = fgetl(fid);
% %             val = str2double(tline);
% %             if isnan(val) || val == -1, break; end
% %             tour(end+1) = val; %#ok<AGROW>
% %         end
% %     end
% % end
% % tour(end+1) = tour(1); % close tour
% % fclose(fid);
% % end
% 
% 
% function [bestTour, bestCost] = tsp_ACO(filename, optTourFile, optVal)
% % Ant Colony Optimization (ACO) for TSP with Candidate Lists
% % filename    : TSPLIB .tsp file
% % optTourFile : (optional) .opt.tour file
% % optVal      : (optional) known optimum value (if no tour file)
% %
% % bestTour    : best tour found
% % bestCost    : cost of best tour
% 
% %% Parameters
% nAnts  = 30;    % number of ants
% nIter  = 1000;   % iterations
% alpha  = 1;     % pheromone influence
% beta   = 5;     % heuristic influence
% rho    = 0.5;   % evaporation rate
% Q      = 100;   % pheromone deposit factor
% candL  = 20;    % candidate list size (neighbors per city)
% 
% %% 1. Read coordinates
% coords = readTSP(filename);
% n = size(coords,1);
% 
% %% 2. Build distance matrix
% D = zeros(n);
% for i = 1:n
%     for j = i+1:n
%         dx = coords(i,1) - coords(j,1);
%         dy = coords(i,2) - coords(j,2);
%         d = sqrt(dx^2 + dy^2);
%         D(i,j) = floor(d+0.5); 
%         D(j,i) = D(i,j);
%     end
% end
% 
% %% 3. Candidate lists (nearest neighbors)
% [~, sortedIdx] = sort(D,2);   % sort neighbors by distance
% candList = sortedIdx(:,2:candL+1); % skip self, take candL nearest
% 
% %% 4. Initialize pheromones
% pheromone = ones(n);
% 
% %% 5. ACO main loop
% bestCost = inf;
% bestTour = [];
% 
% fprintf('Starting ACO on %s with %d ants and %d iterations...\n', filename, nAnts, nIter);
% for it = 1:nIter
%     allTours = zeros(nAnts, n+1);
%     allCosts = zeros(1, nAnts);
% 
%     % construct tours
%     for k = 1:nAnts
%         allTours(k,:) = constructTour(n,D,pheromone,alpha,beta,candList);
%         allCosts(k)   = tourLength(allTours(k,:),D);
%     end
% 
%     % find best in this iteration
%     [minCost, idx] = min(allCosts);
%     if minCost < bestCost
%         bestCost = minCost;
%         bestTour = allTours(idx,:);
%     end
% 
%     % pheromone update
%     pheromone = (1-rho)*pheromone;
%     for k = 1:nAnts
%         Lk = allCosts(k);
%         tour = allTours(k,:);
%         for i = 1:n
%             u = tour(i); v = tour(i+1);
%             pheromone(u,v) = pheromone(u,v) + Q/Lk;
%             pheromone(v,u) = pheromone(u,v);
%         end
%     end
% 
%     % log iteration
%     fprintf('Iter %4d | Best Cost = %d\n', it, bestCost);
% end
% 
% %% 6. Handle optimal value / tour file
% if exist('optTourFile','var') && ~isempty(optTourFile) && isfile(optTourFile)
%     optTour = readOptTour(optTourFile);
%     optCost = tourLength(optTour,D);
%     fprintf('\nACO best cost: %d | Known OPT from tour file: %d\n', bestCost, optCost);
% elseif exist('optVal','var') && ~isempty(optVal)
%     fprintf('\nACO best cost: %d | Provided OPT value: %d\n', bestCost, optVal);
% else
%     fprintf('\nACO best cost: %d (no OPT reference available)\n', bestCost);
% end
% 
% end
% 
% %% ---------------- Helper Functions ----------------
% 
% function coords = readTSP(filename)
% fid = fopen(filename);
% if fid < 0, error('Cannot open file %s',filename); end
% coords = [];
% while true
%     tline = fgetl(fid);
%     if ~ischar(tline), break; end
%     if contains(tline,'NODE_COORD_SECTION'), break; end
% end
% while true
%     tline = fgetl(fid);
%     if ~ischar(tline) || contains(tline,'EOF'), break; end
%     nums = sscanf(tline,'%d %f %f');
%     coords(end+1,:) = nums(2:3)'; %#ok<AGROW>
% end
% fclose(fid);
% end
% 
% function L = tourLength(tour,D)
% L = 0;
% for i = 1:length(tour)-1
%     L = L + D(tour(i),tour(i+1));
% end
% end
% 
% function tour = constructTour(n,D,pheromone,alpha,beta,candList)
% tour = zeros(1,n+1);
% start = randi(n);
% tour(1) = start;
% visited = false(1,n);
% visited(start) = true;
% 
% for i = 2:n
%     u = tour(i-1);
%     candidates = candList(u,:); % nearest neighbors
%     probs = zeros(1,n);
% 
%     % explore candidate list
%     for v = candidates
%         if ~visited(v)
%             tau = pheromone(u,v)^alpha;
%             eta = (1/D(u,v))^beta;
%             probs(v) = tau * eta;
%         end
%     end
% 
%     % fallback: if all candidates are visited, pick from remaining
%     if sum(probs) == 0
%         unvisited = find(~visited);
%         v = unvisited(randi(numel(unvisited)));
%     else
%         probs = probs / sum(probs);
%         v = rouletteWheel(probs);
%     end
% 
%     tour(i) = v;
%     visited(v) = true;
% end
% 
% tour(end) = start; % close the loop
% end
% 
% function v = rouletteWheel(probs)
% r = rand;
% c = cumsum(probs);
% v = find(r <= c,1);
% if isempty(v)
%     v = find(probs>0,1,'last');
% end
% end
% 
% function tour = readOptTour(filename)
% fid = fopen(filename);
% if fid < 0, error('Cannot open file %s',filename); end
% tour = [];
% while true
%     tline = fgetl(fid);
%     if ~ischar(tline), break; end
%     if contains(tline,'TOUR_SECTION')
%         while true
%             tline = fgetl(fid);
%             val = str2double(tline);
%             if isnan(val) || val == -1, break; end
%             tour(end+1) = val; %#ok<AGROW>
%         end
%     end
% end
% tour(end+1) = tour(1); % close tour
% fclose(fid);
% end
% 

function [bestTour, bestCost] = tsp_ACO(filename, optTourFile, optVal)
% Ant Colony Optimization (ACO) for TSP (classical version, no candidate lists)
% filename    : TSPLIB .tsp file
% optTourFile : (optional) .opt.tour file
% optVal      : (optional) known optimum value (if no tour file)
%
% bestTour    : best tour found
% bestCost    : cost of best tour

%% Parameters
nAnts  = 30;    % number of ants
nIter  = 200;   % iterations
alpha  = 1;     % pheromone influence
beta   = 5;     % heuristic influence
rho    = 0.4;   % evaporation rate
Q      = 100;   % pheromone deposit factor

%% 1. Read coordinates
coords = readTSP(filename);
n = size(coords,1);

%% 2. Build distance matrix
D = zeros(n);
for i = 1:n
    for j = i+1:n
        dx = coords(i,1) - coords(j,1);
        dy = coords(i,2) - coords(j,2);
        d = sqrt(dx^2 + dy^2);
        D(i,j) = floor(d+0.5); 
        D(j,i) = D(i,j);
    end
end

%% 3. Initialize pheromones
pheromone = ones(n);

%% 4. ACO main loop
bestCost = inf;
bestTour = [];

fprintf('Starting NORMAL ACO on %s with %d ants and %d iterations...\n', filename, nAnts, nIter);
for it = 1:nIter
    allTours = zeros(nAnts, n+1);
    allCosts = zeros(1, nAnts);

    % construct tours
    for k = 1:nAnts
        tour = constructTour_full(n,D,pheromone,alpha,beta);
        % Optional: apply local 2-opt improvement
        % tour = twoOpt(tour,D);
        allTours(k,:) = tour;
        allCosts(k)   = tourLength(tour,D);
    end

    % find best in this iteration
    [minCost, idx] = min(allCosts);
    if minCost < bestCost
        bestCost = minCost;
        bestTour = allTours(idx,:);
    end

    % pheromone update
    pheromone = (1-rho)*pheromone;
    for k = 1:nAnts
        Lk = allCosts(k);
        tour = allTours(k,:);
        for i = 1:n
            u = tour(i); v = tour(i+1);
            pheromone(u,v) = pheromone(u,v) + Q/Lk;
            pheromone(v,u) = pheromone(u,v);
        end
    end

    % log iteration
    fprintf('Iter %4d | Best Cost = %d\n', it, bestCost);
end

%% 5. Handle optimal value / tour file
if exist('optTourFile','var') && ~isempty(optTourFile) && isfile(optTourFile)
    optTour = readOptTour(optTourFile);
    optCost = tourLength(optTour,D);
    fprintf('\nACO best cost: %d | Known OPT from tour file: %d\n', bestCost, optCost);
elseif exist('optVal','var') && ~isempty(optVal)
    fprintf('\nACO best cost: %d | Provided OPT value: %d\n', bestCost, optVal);
else
    fprintf('\nACO best cost: %d (no OPT reference available)\n', bestCost);
end

end

%% ---------------- Helper Functions ----------------

function coords = readTSP(filename)
fid = fopen(filename);
if fid < 0, error('Cannot open file %s',filename); end
coords = [];
while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    if contains(tline,'NODE_COORD_SECTION'), break; end
end
while true
    tline = fgetl(fid);
    if ~ischar(tline) || contains(tline,'EOF'), break; end
    nums = sscanf(tline,'%d %f %f');
    coords(end+1,:) = nums(2:3)'; %#ok<AGROW>
end
fclose(fid);
end

function L = tourLength(tour,D)
L = 0;
for i = 1:length(tour)-1
    L = L + D(tour(i),tour(i+1));
end
end

function tour = constructTour_full(n,D,pheromone,alpha,beta)
tour = zeros(1,n+1);
start = randi(n);
tour(1) = start;
visited = false(1,n);
visited(start) = true;

for i = 2:n
    u = tour(i-1);
    probs = zeros(1,n);

    for v = 1:n
        if ~visited(v)
            tau = pheromone(u,v)^alpha;
            eta = (1/D(u,v))^beta;
            probs(v) = tau * eta;
        end
    end

    if sum(probs) == 0
        unvisited = find(~visited);
        v = unvisited(randi(numel(unvisited)));
    else
        probs = probs / sum(probs);
        v = rouletteWheel(probs);
    end

    tour(i) = v;
    visited(v) = true;
end

tour(end) = start;
end

function v = rouletteWheel(probs)
r = rand;
c = cumsum(probs);
v = find(r <= c,1);
if isempty(v)
    v = find(probs>0,1,'last');
end
end

function tour = twoOpt(tour,D)
improved = true;
while improved
    improved = false;
    n = length(tour)-1; % last = start
    for i = 2:n-2
        for j = i+1:n-1
            d1 = D(tour(i-1),tour(i)) + D(tour(j),tour(j+1));
            d2 = D(tour(i-1),tour(j)) + D(tour(i),tour(j+1));
            if d2 < d1
                tour(i:j) = tour(j:-1:i);
                improved = true;
            end
        end
    end
end
end

function tour = readOptTour(filename)
fid = fopen(filename);
if fid < 0, error('Cannot open file %s',filename); end
tour = [];
while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    if contains(tline,'TOUR_SECTION')
        while true
            tline = fgetl(fid);
            val = str2double(tline);
            if isnan(val) || val == -1, break; end
            tour(end+1) = val; %#ok<AGROW>
        end
    end
end
tour(end+1) = tour(1);
fclose(fid);
end
