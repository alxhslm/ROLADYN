function [Qi,Ki] = race_compliance_loads(Race,wi)
Qi = (Race.Kax * wi / Race.R + Race.Kfl * wi ./ (Race.R + wi).^3);
Ki = Race.Kax/Race.R + Race.Kfl*(1./(Race.R + wi).^3 + 3*wi./(Race.R + wi).^4);