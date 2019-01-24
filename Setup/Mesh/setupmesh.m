function P = setupmesh(P)
NDofe = 4;

%find the number of rotor dof
NRotor = length(P.Rotor);
iDofCount = 0;
for i = 1:NRotor
    P.Rotor{i}.Bearing = {};
    NDofRotor(i) = P.Rotor{i}.NDof;
    P.Rotor{i}.iGlobal = iDofCount + (1:P.Rotor{i}.NDof);
    iDofCount = iDofCount + P.Rotor{i}.NDof;
end
P.Mesh.Rotor.NDof = sum(NDofRotor);

%and the bearings
NBearings = length(P.Bearing);
for i = 1:NBearings
    P.Bearing{i}.iGlobal = iDofCount + (1:NDofe);
    iDofCount = iDofCount + NDofe;
end
P.Mesh.Bearing.NDof = NDofe * NBearings;

%work out the number of internal states
NInternal = zeros(NBearings,2);
for i = 1:NBearings
    NInternal(i,:) = P.Bearing{i}.NDofInt;
end
P.Mesh.Bearing.NDofInt = sum(NInternal(:));
P.Mesh.NDofInt = P.Mesh.Bearing.NDofInt;
 
for i = 1:length(P.Excite)
    NExcite(i) = P.Excite{i}.NInput;
end
P.Mesh.Excite.NInput = sum(NExcite(:));

P.Mesh.NInput = P.Mesh.Excite.NInput;
P.Mesh.NDof = P.Mesh.Rotor.NDof + P.Mesh.Bearing.NDof;
P.Mesh.NDofTot = P.Mesh.Rotor.NDof + P.Mesh.Bearing.NDof + P.Mesh.NDofInt;
