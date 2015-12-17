%function to make putting together the string less of a pain
function [content] = AddLine(local_s,s)
    local_s = [local_s sprintf('%s \n',s)];    
    content = local_s;
end