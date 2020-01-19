function [sigma,z] = shaft_stress(P,X)
for i = 1:length(P.Rotor)
   Sr = P.Rotor{i}.S;
   for j = 1:length(P.Rotor{i}.Shaft)
       Ss = P.Rotor{i}.Shaft{j}.S * Sr;
       Z = P.Rotor{i}.Shaft{j}.Section.Z;
       
       for k = 1:length(P.Rotor{i}.Shaft{j}.Element)
           Se = P.Rotor{i}.Shaft{j}.Element{k}.S*Ss;
           Re = P.Rotor{i}.Shaft{j}.Element{k}.R;
            
           Fs = mtimesx(Re' * P.Rotor{i}.Shaft{j}.Element{k}.K * Re * Se,X);
           My = Fs(3,:,:);
           Mx = Fs(4,:,:);
       
           M = hypot(Mx,My);
           for n = 1:size(M,2)
                sigma{i}{j}(k,n,:) = M(:,n,:) / Z;
           end
           z{i}{j}(k) = P.Rotor{i}.Shaft{j}.Element{k}.z;
       end
   end
end

