function [roi] = parimgbox_mix(image, bbox)
    roi = image(bbox(2) : bbox(2) + bbox(4), bbox(1) : bbox(1) + bbox(3), :);
end