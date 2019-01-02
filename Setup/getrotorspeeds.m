function N = getrotorspeeds(R)
for i = 1:length(R)
    N(i) = R{i}.Speed;
end
