function [bestTour, bestCost] = tsp_DP_bitmask(filename)
% Exact TSP solution using Dynamic Programming + Bitmasking
% filename : TSPLIB .tsp file
% bestTour : optimal tour
% bestCost : optimal cost

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

%% 3. DP table: dp(mask, i) = min cost to visit nodes in mask ending at i
INF = 1e9;
dp = INF * ones(2^n, n);
parent = zeros(2^n, n);

% starting at node 1
dp(1,1) = 0;

%% 4. Fill DP table
for mask = 1:2^n-1
    for u = 1:n
        if bitand(mask, 2^(u-1)) == 0
            continue; % u not in mask
        end
        for v = 1:n
            if bitand(mask, 2^(v-1)) ~= 0
                continue; % v already visited
            end
            nextMask = bitor(mask, 2^(v-1));
            cost = dp(mask,u) + D(u,v);
            if cost < dp(nextMask,v)
                dp(nextMask,v) = cost;
                parent(nextMask,v) = u;
            end
        end
    end
end

%% 5. Close tour and find min cost
fullMask = 2^n - 1;
bestCost = INF;
lastNode = 1;
for i = 2:n
    cost = dp(fullMask,i) + D(i,1);
    if cost < bestCost
        bestCost = cost;
        lastNode = i;
    end
end

%% 6. Reconstruct tour
bestTour = zeros(1,n+1);
mask = fullMask;
node = lastNode;
for i = n:-1:1
    bestTour(i) = node;
    prev = parent(mask,node);
    mask = bitxor(mask, 2^(node-1));
    node = prev;
end
bestTour(end) = bestTour(1); % close tour

fprintf('DP TSP optimal cost: %d\n', bestCost);

end

%% ---------------- Helper Function ----------------
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
