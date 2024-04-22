function [M] = get_cleanMask(I,M,maskFile)
    M = imcomplement(M);
    M = getWatershedMask(M,false,6,12);
    M = logical(M);
    M = bwareafilt(M,[30 300]);
    M1 = getWatershedMask(I,true,6,12);
    M1 = logical(M1);
    M1 = bwareafilt(M1,[30 300]);
    M = M+M1;
    imwrite(M,maskFile);
end