% This is 

function [vector] = matrix2vector (matrix, condition)

cnt = 1;

if condition == "full"
    
    for i = 1:size(matrix,1)
        for j = 1:size(matrix,2)
            vector(cnt) = matrix(i,j);
            cnt = cnt + 1;
        end
    end
        
elseif condition == "lower"
    
    for i = 2:size(matrix,1)
        for j = 1: (i-1)
            vector(cnt) = matrix(i,j);
            cnt = cnt + 1;
        end
    end
    
end

end