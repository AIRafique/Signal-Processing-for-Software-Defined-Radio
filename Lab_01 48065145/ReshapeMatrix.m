function [outputMatrix] = ReshapeMatrix(size)
    array = 0:size^2 - 1;
    outputMatrix = reshape(array, [size, size]);
end