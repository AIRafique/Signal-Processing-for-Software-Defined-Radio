function [OutputMatrix] = ComputationMatrix(size)
    OutputMatrix = zeros(size,size);
    arrayA = 0:size-1;
    for i=1:size
        OutputMatrix(i,:) = arrayA;
        arrayA = arrayA+size;
    end
    OutputMatrix = OutputMatrix';
end