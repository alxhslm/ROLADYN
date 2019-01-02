function F = packharm(Fh,NHarm)
F = zeros(NHarm+1,length(Fh));          
F(2,:) = Fh(:);