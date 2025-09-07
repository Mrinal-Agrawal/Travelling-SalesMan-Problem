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
