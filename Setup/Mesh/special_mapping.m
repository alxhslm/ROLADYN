function S = special_mapping(NDofe,NDofTot,iGlobal,iNode)
S = zeros(NDofe,NDofTot);
S((1:NDofe),(iGlobal(iNode) - 1)*NDofe+(1:NDofe)) = eye(NDofe);