function output = elevcut(input,elevation,elevation_cutoff)
% elevation angle cutoff
% Inputs:
%          input           = matrix input (C1, Sat_position)
%          elevation       = elevation angle
%          elvation_cutoff = elevation mask (set default 15)
% Outputs:
%          output          = matrix output, that already cut off.

elevation = abs(elevation);
mask                         = elevation;
mask(mask<elevation_cutoff)  = NaN;
mask(~isnan(mask))           = 1;
    
[s1, s2] = size(input);
    
for i = 1:s1
    mask_s(i,:) = mask;
end
    
inputs = mask_s.*ones(s1,s2).*input;

for cn = 1:s1         % Clear nan value
    cnn = inputs(cn,:);
    cnn = cnn(~isnan(cnn));
    output(cn,:)=cnn;  
end

end