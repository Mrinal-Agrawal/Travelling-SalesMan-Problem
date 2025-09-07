function cost = evaluateTour(tour, D)
n = length(tour);
cost = 0;
for i=1:n-1
    cost = cost + D(tour(i), tour(i+1));
end
cost = cost + D(tour(end), tour(1));
end
