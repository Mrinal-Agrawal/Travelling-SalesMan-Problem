function [tour, tourCost] = tsp_MST_approx(filename)
%TSP_MST_APPROX 2-approximation TSP using MST preorder DFS
coords = readTSP(filename);
n = size(coords,1);

%% 1. Build distance matrix (TSPLIB EUC_2D rounding)
D = zeros(n);
for i = 1:n
    for j = i+1:n
        dx = coords(i,1) - coords(j,1);
        dy = coords(i,2) - coords(j,2);
        d = sqrt(dx^2 + dy^2);
        D(i,j) = floor(d + 0.5); % TSPLIB rounding
        D(j,i) = D(i,j);
    end
end

%% 2. Compute MST
mstAdj = primMST(D);

%% 3. DFS preorder traversal of MST
visited = false(1,n);
dfsTour = [];
dfsPreorder(1); % start DFS at node 1

%% 4. Shortcut repeated vertices
[tour, tourCost] = shortcutTour(dfsTour, D);

%% 5. Plot the tour
figure;
plot(coords(:,1), coords(:,2), 'ro', 'MarkerSize',8,'MarkerFaceColor','r'); hold on;
for i=1:n
    text(coords(i,1)+10, coords(i,2)+10, num2str(i));
end
tourCoords = coords(tour,:);
plot(tourCoords(:,1), tourCoords(:,2), 'b-', 'LineWidth',1.5);
title(sprintf('2x-MST Approx TSP Tour, Cost = %.0f', tourCost));
axis equal; grid on;

%% ---------------- Nested helper ----------------
    function dfsPreorder(u)
        visited(u) = true;
        dfsTour(end+1) = u; %#ok<AGROW>
        for v = mstAdj{u}
            if ~visited(v)
                dfsPreorder(v);
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
    if contains(tline,'NODE_COORD_SECTION')
        break;
    end
end
while true
    tline = fgetl(fid);
    if ~ischar(tline) || contains(tline,'EOF'), break; end
    nums = sscanf(tline,'%d %f %f');
    coords(end+1,:) = nums(2:3)'; %#ok<AGROW>
end
fclose(fid);
end

function mstAdj = primMST(D)
n = size(D,1);
inMST = false(1,n);
key = inf(1,n);
parent = zeros(1,n);
key(1) = 0;
for cnt=1:n
    temp = key; temp(inMST) = inf;
    [~, u] = min(temp);
    inMST(u) = true;
    for v=1:n
        if ~inMST(v) && D(u,v) < key(v)
            key(v) = D(u,v);
            parent(v) = u;
        end
    end
end
mstAdj = cell(n,1);
for v=2:n
    u = parent(v);
    mstAdj{u}(end+1) = v;
    mstAdj{v}(end+1) = u;
end
end

function [tour, cost] = shortcutTour(nodeSeq, D)
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
% Close the tour
cost = cost + D(tour(end), tour(1));
tour(end+1) = tour(1);
end
