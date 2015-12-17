function [varnames2,varvalue2,vartypes,count] = Readsettings(name,nrproc)        
        %mkdir(['dir_' name]);
        %cd (['dir_' name]);
        %cmd = ['"C:\Program Files\7-Zip\7z.exe" e ../' name '.tar.gz'];        
        %[status,result1] = system(cmd );
        %cmd = ['"C:\Program Files\7-Zip\7z.exe" e ' name '.tar'];
        cmd = ['tar -xf ' name '.tar.gz'];        
        [status,result2] = system(cmd );
        cd (name);
        if (nrproc > 1)
            fid = fopen('out_settings_0.txt','r');
        else
            fid = fopen('out_settings.txt','r');
        
        end
        varstat = 1;
        count = 0;
        varnames = {};
        vartypes = {};
        varvalues = {};
        tline = fgetl(fid);
        while ischar(tline)
            t_line =  strtrim(tline);
            if (varstat == 1)
                if (~isempty(t_line))
                    count = count + 1;
                    varnames{count} = strtrim(tline);
                    varstat = 0;
                end
            else
                value = t_line;
                x = str2num(value);
                if (isempty(x))
                    if (strcmpi(value,'F') || strcmpi(value,'T'))
                        vartypes{count} = 'int(1)';
                        if (strcmpi(value,'F'))
                            varvalue{count} = '0';
                        else
                            varvalue{count} = '1';
                        end
                    else
                        vartypes{count} = 'varchar(300)';
                        varvalue{count} = value;
                    end
                else
                    vartypes{count} = 'double';                   
                    varvalue{count} = num2str(str2num(value));
                    if ( strcmpi( varvalue{count},'-Inf') )
                        varvalue{count} = '0';
                    end
                    if ( strcmpi( varvalue{count}, 'Inf') )
                        varvalue{count} = '0';
                    end
                end
                
                varstat = 1;
            end
            tline = fgetl(fid);
        end
        fclose(fid);
        varnames2 = strtrim(varnames);
        varvalue2 = strtrim(varvalue);
end