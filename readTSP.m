function [coords, distMatrix] = readTSP(filename)
%READTSP Reads a TSP file and returns coordinates and distance matrix
%   Supports EDGE_WEIGHT_TYPE = EUC_2D and ATT

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    nCities = 0;
    coords = [];
    EDGE_WEIGHT_TYPE = 'EUC_2D'; % default

    % --- Read header and coordinates ---
    while true
        line = fgetl(fid);
        if line == -1
            error('Unexpected end of file.');
        end

        if contains(line, 'DIMENSION')
            tokens = regexp(line,'\d+','match');
            nCities = str2double(tokens{1});
        elseif contains(line, 'EDGE_WEIGHT_TYPE')
            parts = strsplit(line, ':');
            EDGE_WEIGHT_TYPE = strtrim(parts{2});  % Correctly extract type
        elseif contains(line, 'NODE_COORD_SECTION')
            coords = zeros(nCities,2);
            for i = 1:nCities
                data = fscanf(fid, '%d %f %f', 3);
                coords(i,:) = data(2:3);
            end
            break;
        end
    end
    fclose(fid);

    % --- Compute distance matrix ---
    distMatrix = zeros(nCities);
    for i = 1:nCities
        for j = 1:nCities
            if i ~= j
                dx = coords(i,1) - coords(j,1);
                dy = coords(i,2) - coords(j,2);

                switch EDGE_WEIGHT_TYPE
                    case 'EUC_2D'
                        distMatrix(i,j) = round(sqrt(dx^2 + dy^2));
                    case 'ATT'
                        rij = sqrt((dx^2 + dy^2)/10.0);
                        tij = round(rij);
                        if tij < rij
                            distMatrix(i,j) = tij + 1;
                        else
                            distMatrix(i,j) = tij;
                        end
                    otherwise
                        error('Unsupported EDGE_WEIGHT_TYPE: %s', EDGE_WEIGHT_TYPE);
                end
            end
        end
    end
end

