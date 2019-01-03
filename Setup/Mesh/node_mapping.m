function S = node_mapping(NDofe,NDofTot,iGlobal)
S = zeros(NDofe,NDofTot);
S((1:NDofe),iGlobal+(1:NDofe)) = eye(NDofe);
