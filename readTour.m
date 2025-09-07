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
