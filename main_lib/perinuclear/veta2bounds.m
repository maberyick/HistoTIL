function bounds = veta2bounds(ep_nuclei)
%parpool(4)
% transform into bounds struct
bounds = [];
for i = 1:numel(ep_nuclei)
    x = ep_nuclei{i};
    bounds(i).r = x(:,1);
    bounds(i).c = x(:,2);
    [bounds(i).centroid_r,bounds(i).centroid_c] = poly_centroid(x(:,1),x(:,2));
end

% check bounds struct
for i = 1:numel(bounds)
    if size(bounds(end).r,2) == 1
        bounds(i).r = bounds(i).r';
        bounds(i).c = bounds(i).c';
    end
end