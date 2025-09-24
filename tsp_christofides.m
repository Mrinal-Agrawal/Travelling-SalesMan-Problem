function [tour, tourCost] = tsp_christofides(filename, varargin)
% Christofides algorithm for metric TSP (EUC_2D)
% filename : TSPLIB .tsp file
% varargin : optional: either optTourFile or numeric optCost
% tour     : Christofides tour node sequence
% tourCost : Christofides tour total cost

%% 1. Read coordinates
coords = readTSP(filename);
n = size(coords,1);

%% 2. Build distance matrix (TSPLIB EUC_2D rounding)
D = zeros(n);
for i = 1:n
    for j = i+1:n
        dx = coords(i,1) - coords(j,1);
        dy = coords(i,2) - coords(j,2);
        d = sqrt(dx^2 + dy^2);
        D(i,j) = floor(d + 0.5);
        D(j,i) = D(i,j);
    end
end

%% 3. Compute MST via Prim's algorithm
mstAdj = primMST(D);

%% 4. Find nodes with odd degree
deg = cellfun(@numel, mstAdj);
oddNodes = find(mod(deg,2)==1);

%% 5. Minimum-weight perfect matching of odd-degree nodes
oddD = D(oddNodes, oddNodes);
nOdd = length(oddNodes);
for i = 1:nOdd
    oddD(i,i) = 1e6;  % large finite number to prevent self-match
end
[matches, ~] = matchpairs(oddD, 1e6);

% Map back to original node indices
matching = zeros(size(matches));
for k = 1:size(matches,1)
    matching(k,:) = oddNodes(matches(k,:));
end

%% 6. Form multigraph (MST + matching)
multiAdj = mstAdj;
for k = 1:size(matching,1)
    u = matching(k,1); v = matching(k,2);
    multiAdj{u}(end+1) = v;
    multiAdj{v}(end+1) = u;
end

%% 7. Eulerian tour (DFS)
% visited = false(1,n);
eulerTour = [];
dfsEuler(1);

%% 8. Shortcut repeated vertices to get TSP tour
[tour, tourCost] = shortcutTour(eulerTour, D);

%% 9. Process optional input (optTourFile or optCost)
optCost = [];
if ~isempty(varargin)
    arg = varargin{1};
    if ischar(arg) % it's a .opt.tour file
        optTour = readOptTour(arg);
        optCost = 0;
        for i = 1:length(optTour)-1
            optCost = optCost + D(optTour(i), optTour(i+1));
        end
    elseif isnumeric(arg) % it's the optimal cost
        optCost = arg;
    end
end

%% 10. Plot Christofides tour
figure; hold on; grid on; axis equal;
plot(coords(:,1), coords(:,2), 'ko','MarkerSize',6,'MarkerFaceColor','k');
tourCoords = coords(tour,:);
plot(tourCoords(:,1), tourCoords(:,2), 'b-', 'LineWidth',1.5, ...
    'DisplayName', sprintf('Christofides: %.0f', tourCost));

% Optionally show optimal tour or cost in legend
if ~isempty(optCost)
    if exist('optTour','var')
        optCoords = coords(optTour,:);
        plot(optCoords(:,1), optCoords(:,2), 'r--', 'LineWidth',1.5, ...
            'DisplayName', sprintf('Optimal: %.0f', optCost));
    else
        % Only show cost in legend
        plot(NaN, NaN, 'r--', 'DisplayName', sprintf('Optimal cost: %.0f', optCost));
    end
end

legend('show');
title('TSP Tour Comparison');

fprintf('Christofides tour cost: %.0f\n', tourCost);
if ~isempty(optCost)
    fprintf('Optimal tour cost: %.0f\n', optCost);
end

%% Nested DFS for Eulerian tour
    function dfsEuler(u)
        eulerTour(end+1) = u;
        for v = multiAdj{u}
            if any(multiAdj{v}==u)
                multiAdj{v}(multiAdj{v}==u) = [];
                multiAdj{u}(multiAdj{u}==v) = [];
                dfsEuler(v);
            end
        end
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

function tour = readOptTour(filename)
fid = fopen(filename);
if fid < 0, error('Cannot open file %s',filename); end
tour = [];
while true
    tline = fgetl(fid);
    if ~ischar(tline), break; end
    if contains(tline,'TOUR_SECTION'), break; end
end
while true
    tline = fgetl(fid);
    if ~ischar(tline) || contains(tline,'-1') || contains(tline,'EOF'), break; end
    val = sscanf(tline,'%d');
    tour(end+1) = val; %#ok<AGROW>
end
tour(end+1) = tour(1); % close tour
fclose(fid);
end

function mstAdj = primMST(D)
n = size(D,1);
inMST = false(1,n);
key = inf(1,n);
parent = zeros(1,n);
key(1) = 0;
for cnt = 1:n
    temp = key; temp(inMST) = Inf;
    [~, u] = min(temp);
    inMST(u) = true;
    for v = 1:n
        if ~inMST(v) && D(u,v) < key(v)
            key(v) = D(u,v);
            parent(v) = u;
        end
    end
end
mstAdj = cell(n,1);
for v = 2:n
    u = parent(v);
    mstAdj{u}(end+1) = v;
    mstAdj{v}(end+1) = u;
end
end

function [tour, cost] = shortcutTour(nodeSeq,D)
visited = false(1,size(D,1));
tour = [];
cost = 0;
prev = -1;
for u = nodeSeq
    if ~visited(u)
        tour(end+1) = u; %#ok<AGROW>
        visited(u) = true;
        if prev ~= -1
            cost = cost + D(prev,u);
        end
        prev = u;
    end
end
cost = cost + D(tour(end),tour(1));
tour(end+1) = tour(1);
end
