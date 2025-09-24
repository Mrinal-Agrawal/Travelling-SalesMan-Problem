function optCost = tsp_optTourCost(tspFile, tourFile)
    % Read coordinates
    coords = readTSP(tspFile);
    n = size(coords,1);

    % Build distance matrix
    D = zeros(n);
    for i=1:n
        for j=i+1:n
            dx = coords(i,1)-coords(j,1);
            dy = coords(i,2)-coords(j,2);
            d = sqrt(dx^2 + dy^2);
            D(i,j) = floor(d+0.5); % TSPLIB rounding
            D(j,i) = D(i,j);
        end
    end

    % Read tour file
    fid = fopen(tourFile);
    tour = [];
    while true
        tline = fgetl(fid);
        if ~ischar(tline), break; end
        if contains(tline,'TOUR_SECTION'), break; end
    end
    while true
        tline = fgetl(fid);
        if ~ischar(tline) || contains(tline,'-1') || contains(tline,'EOF'), break; end
        tour(end+1) = str2double(tline); %#ok<AGROW>
    end
    fclose(fid);

    % Compute tour cost
    optCost = 0;
    for i=1:length(tour)-1
        optCost = optCost + D(tour(i), tour(i+1));
    end
    optCost = optCost + D(tour(end), tour(1)); % close the tour
end

function coords = readTSP(filename)
    fid = fopen(filename);
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
