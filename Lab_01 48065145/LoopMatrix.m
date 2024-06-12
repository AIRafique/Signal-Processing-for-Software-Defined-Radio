function [outputMatrix] = LoopMatrix(size)
outputMatrix = zeros(size, size);
    for i = 1:size
        for j=1:size
            outputMatrix(i,j) = (i-1) + (j-1)*size;
        end
    end
end