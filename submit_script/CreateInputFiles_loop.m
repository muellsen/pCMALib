function [vars,outputs,out_size] =  CreateInputFiles(opts,defopts)
%% check which parameters are getting varied
namec = fieldnames(opts);
pers = ones(size(namec,1),1);

for i = 1:size(namec,1)
    name = namec{i};    
    if (~(strcmp(name,'CECFolders') || strcmp(name,'seed_folder') || strcmp(name,'output_folder') || strcmp(name,'funcName')))
        if (size(myeval(getfield(opts,name)),2) > 1)
            pers(i) = size(myeval(getfield(opts,name)),2);            
        end
    end    
end

out_size = prod(pers);
outputs = cell(size(namec,1),out_size);
sw = zeros(out_size,1);
vars = cell(size(namec,1),1);
for i = 1:size(namec,1)
    per_c = 1;
    c_string = 0;
    boolean = 0;
    name = namec{i}; 
    vars{i} = upper(name);
    if (~(strcmp(name,'CECFolders') || strcmp(name,'seed_folder') || strcmp(name,'output_folder') || strcmp(name,'funcName')))
         value = myeval(getfield(opts,name));
         defvalue = getfield(defopts,name);
         if ( ~(isempty(findstr('true',defvalue)) && isempty (findstr('false',defvalue))))
            boolean = 1; 
         end
    else
         value = getfield(opts,name);
         pos = strfind(value, '%');
        if (~isempty(pos))
            value = value(1:pos(1)-1);
        end
         c_string = 1;
    end       
    if (pers(i) > 1)
       last_div = find(sw,1); 
       if (isempty(last_div))
           last_div = out_size;
       end
       x = last_div/pers(i);
       sw(x) = 1;
    end
    count = 1;
    for (j = 1:out_size)        
        if (c_string)
            outputs{i,j} = value;                 
        else
            %we have to replace 0/1 with true false for booleans
            if (boolean)
                outputs{i,j} = value(per_c);
                if (outputs{i,j})
                    outputs{i,j} = 'TRUE';
                else
                    outputs{i,j} = 'FALSE';
                end
            else
                outputs{i,j} = num2str(value(per_c));
            end
        end
        count = count + 1;
        if (count > find(sw,1)) 
           per_c = per_c + 1;
           if (per_c > pers(i))
               per_c = 1;
           end
        count = 1;
        end
        
    end
end
end