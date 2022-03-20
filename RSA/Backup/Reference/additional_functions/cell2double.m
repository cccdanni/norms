function [x] = cell2double(cellarray)

x=zeros(size(cellarray,1),size(cellarray,2));

for z=1:size(cellarray,1)
    for s=1:size(cellarray,2)
        x(z,s)=str2double(cell2mat(cellarray(z,s)));
    end
end