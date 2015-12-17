%collects multiple Runs together and submits them once they fill a certain
%time window
function [nr_o] = AddtoSub(pathtoCMA,queuetime,expt,out,mpi)
persistent timesize;
persistent l_count;
persistent batch;
persistent nr;

%initialize
if ((strcmpi(pathtoCMA,'')))
   timesize = 0;
   nr = 0;
   l_count = 0;
   batch = '#!/bin/sh';
   batch = sprintf('%s\n',batch);
   %batch = [batch 'bsub -W 08:00 -o report' num2str(nr) ' << EOF'];
   %batch = sprintf('%s\n',batch);
else
    s = '';
    s = sprintf('%s;\n',s);
    if ( (timesize + expt > queuetime || l_count > 600 || expt == 0) && l_count > 0) %we have a full batch 
                                         %submit it
        fid = fopen(strcat(num2str(nr),'batch.txt'),'wt');
        fprintf(fid,'%s',batch);
        fclose(fid);
        nr = nr + 1;
        timesize = 0;
        l_count = 0;
        batch = '#!/bin/sh';
        batch = sprintf('%s\n',batch);
        %batch = [batch 'bsub -W 08:00 -o report' num2str(nr) ' << EOF'];
        %batch = sprintf('%s\n',batch);
    end
   timesize = timesize + expt;
   l_count = l_count + 6;

   %extract input


   %disp(input_file);
   batch = [batch 'mkdir ' out s];
   batch = [batch 'tar -x ' out '.txt -f inputs.tar.gz -O > ' out '/input.txt' s];         
   if (mpi > 1)
    batch = [batch 'ompirun ' pathtoCMA ' ' out '/input.txt > ' out '/stdout.txt' s];   
   else
    batch = [batch pathtoCMA ' ' out '/input.txt > ' out '/stdout.txt' s];
   end
   batch = [batch 'tar czvf ' out '.tar.gz ' out s];
   batch = [batch 'rm -r ' out s];        

end

nr_o = nr;
end