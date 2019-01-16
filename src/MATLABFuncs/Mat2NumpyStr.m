function [str]=Mat2NumpyStr(mat)  

str=num2str(mat,'%.1f,');
for i=1:size(str,1)
    str(i,end)=';';
end
str = reshape(str.',1,[]);
str=['[[' str ']]'];
str = strrep(str,';','],[');
str = strrep(str,',[]','');
str = strrep(str,' ','');