function [something] = Brutus_Benchmarks_final(pathtoCMA,workdir,name,repeats,inopts,queuetime,indbname,inexpt,innrproc)
%for this to work brutus must be reachable via ssh without password
%(keypair) and the jdbc driver has to be access able (most likly the path
%must be adjusted)



%% DEF Options
defopts.TolFun          = '1e-12 % stop if fun-changes smaller TolFun';
defopts.dimensions      = '10 % the dimension of the problem';
defopts.rel_sigma       = '0.2 % size of the inital sigma relative to the bounds - only used if ABS_SIGMA is not set';
defopts.abs_sigma       = '0 % size of the inital $\sigma$ as absolut value';
defopts.alldim_lBounds  = '-Inf % Lower bounds on all dimensions';
defopts.alldim_uBounds  = 'Inf %Upper bounds on all dimensions';
defopts.PopSize         = '0%(4 + floor(3*log(N)))  population size, lambda';
defopts.use_CEC         = 'false % use the CEC2005 benchmark suite';
defopts.use_LJC         = 'false %use the Lennard Jones potential';
defopts.use_TIP         = 'false % use the TIP4P water potential as target function';
defopts.use_BBOB        = 'false %  use the BBOB benchmark suite as target function';
defopts.Benchfctnr      = '1 %which function of the CEC2005 or BBOB benchmark suite to use';
defopts.output_folder   = 'out % output data is saved into this folder (relative to workdir)';
defopts.StopFitness     = '-Inf %stop if f(xmin) < stopfitness, minimization';
defopts.StopMaxFunEvals = 'Inf %maximal number of fevals';
defopts.StopMaxIter     = '0 % 1e3*(N+5)^2/sqrt(popsize) maximal number of iterations';
defopts.TolX            = '0 %1e-11*max(insigma) % stop if x-change smaller TolX';
defopts.TolUpX          = '0 %1e3*max(insigma) % stop if x-changes larger TolUpX';
defopts.TolFun          = '1e-12 % stop if fun-changes smaller TolFun';
defopts.TolHistFun      = '1e-13 % stop if back fun-changes smaller TolHistFun';
defopts.StopOnWarnings  = 'true  % ''no''==''off''==0, ''on''==''yes''==1 ';
defopts.WarnOnEqualFunctionValues = 'true  % ';
defopts.EvalInitialX    = 'true  % evaluation of initial solution';
defopts.ParentNumber    = '0 %floor(popsize/2)     % popsize equals lambda';
defopts.RecombinationWeights = '3 % Super--linear (3), linear (2) or equal (1)';
defopts.VerboseModulo   = '100 %Messaging after every i-th iteration';
defopts.flgGenData      = 'false % Save data of all iterations';
defopts.intGenData      = '1 %generation interval to log data, default is interval = 1';
defopts.funcName        = 'no name assigned %Objective Function Name';
defopts.CECFolders      = ' %where to find the supportData folder for the CEC2005 benchmark suite relative to the working directory';
defopts.silentMPI       = 'true %if only rank 0 should report output';
defopts.pscma           = 'false % Switch PS-CMA-ES on or off';
defopts.psoWeight       = '0.7 %how much the PSO update influences the covariance';
defopts.psoFreq         = '200 %Intervall length between PSO updates';
defopts.accuracy        = '0 % Successful run if global min - f(x) < accuracy';
defopts.global_min      = '-Inf %Global minimum';
defopts.use_seed        = 'false %if a seed should be used';
defopts.seed_folder     = '%folder containing the seed file seed.txt'; 
defopts.qr_sampling     = 'true %usage of Quasi Random sampling';
defopts.qr_sampler      = '1 %0.............Sobol, (default)1 Sobol with scrambling, 2 Halton,3 Halton R implementation, 4 Faure (buggy!), 5 Niederreiter';
defopts.qr_inverter     = '1 %(0)Moros Inverse\\&&(1)Peter J. Acklam?s Inverter\\&&(2)Inverter from the R';
defopts.use_init_bounds = 'false % if special bounds should be used for initalization of population';
defopts.init_uBounds    = 'Inf %Upper bounds for initialization';
defopts.init_lBounds    = '-Inf % Lower bounds for initialization';
defopts.benchmark       = 'false % Switch Benchmark on or off, this causes CMA to record and keep track of several variables that are required in the benchmark protocols';
defopts.record_accuracy = '0.0 %record when the optimizer reaches this level of accuracy';
defopts.restart_cma     = 'false % if multi restart CMA (IPOP) should be used ';
defopts.restart_type    = '0 % (0) restart randomly within bounds\\&& (1) restart from point of convergence\\&& (2) restart from same startpoint all the time';
defopts.Restarts        = '0 % limit on how many restarts are allowed, 0 = unlimited ';
defopts.IncPopSize      = '1.25 % factor by which the population size is increased every restart ';
defopts.record_besthist = 'false % record a history of the fitness over time ';
defopts.record_modulo   = '100 %!how many recordings (equal spaced) of the bestever should be made';
defopts.MaxIncFac       = '10  %! this setting is important because for big problem sizes one easily gets into memory troubles';
defopts.BFGS_use        = 'false % if BFGS should be used to assist CMA. This is still in development!';
defopts.BFGS_position   = '2 % 1 = replace X values by local minimum X\\&& 2 = replace F values with F values at local minimum';
defopts.BFGS_factr      = '0 %if zero StopTolFun/machine eps is used';
defopts.BFGS_pgtol      = '1.0d-5 %gradient quality stop criteria';
defopts.BFGS_central_difference = 'false %if central difference should be used, otherwise backward difference is utilized';
defopts.BFGS_confprop   = '0.95 %BFGS_confprop is the upper limit, on how far BFGS is allowed to travel within the searchspace relative to the current density of the multivariant gauss of CMA adjusted in terms of probability (e.g. 0.95 means that BFGS is only allowed to travel within the 0.95% probability volumn of the multivariant gauss)';
defopts.BFGS_sigmastep_min = '0.1 %! BFGS_sigmastep_min is the lower limit, on how far BFGS has to travel within each iteration relative to the current sigma';
defopts.StopTimeUse     = 'false % usage of a stop time';
defopts.StopTimeHH      = '0 %hours';
defopts.StopTimeMM      = '0 %minutes';
defopts.StopTimeSS      = '0 %seconds';
defopts.write_pdb       = 'false % if a pdb file representing the current best solution is written in intervals according to the';
defopts.LJ_comp         = '1 % compression term for the Lennard Jones potential with compression';
defopts.use_LJ_comp     = 'false % use the Lennard Jones potential with compression as target function';
defopts.use_MATFUNC     = 'false % use the template for a Matlab target as function';
defopts.flgouttxt       = 'true %saves the results of the pCMAlib run into textfiles';
defopts.flgGenTrace     = 'false %log data from the trace (bestever.f and bestever.x) in intervals given by intGenData, this is included in flgGenData';
defopts.use_DF          = 'false % use the DoubleFunnel benchmark as target function';
defopts.DF_rast         = 'true % if rastrigin function should be applied to the Double Funnel benchmark';
defopts.DF_s            = '0 %setting for parameter s for the Double Funnel benchmark';
defopts.use_RANDOM_LANDSCAPE = 'false %if a random landscape should be used as test function';







%% ---------------------- Handling Input Parameters -----------------------
if nargin < 1 || isequal(pathtoCMA, 'defaults') % pass default options
  if nargin < 1
    disp('Default options returned (type "help Brutus_Benchmarks_final" for help).');
  end
  something = defopts;
  if nargin > 1 % supplement second argument with default options
    something = getoptions(workdir, defopts);
  end
  return;
end
if isequal(pathtoCMA, 'displayoptions')
 names = fieldnames(defopts); 
 for name = names'
   disp([name{:} repmat(' ', 1, 20-length(name{:})) ': ''' defopts.(name{:}) '''']); 
 end
 return; 
end
% Compose options opts
if nargin < 5 || isempty(inopts) % no input options available
  inopts = []; 
  opts = defopts;
else
  opts = getoptions(inopts, defopts);
end

if nargin < 7 || isempty(indbname) % no input options available
  dbname = name;
else
  dbname = indbname;
end
if nargin < 8 || isempty(inexpt) % if time for jobs is given
  expt = 0; %every job its own batch
else
  expt = inexpt;
end
if nargin < 9 || isempty(innrproc) % single job or multiple processors
  nrproc = 1; %single job default
else
  nrproc = innrproc;
end       



%%-----------------Checking the status of the task-------------------------
cmd = ['ssh brutus.ethz.ch ''ls ' workdir name ''''];
[status,result] = system(cmd);
if ( findstr(result,'No such file or directory'))
    cmd = ['ssh brutus.ethz.ch ''mkdir ' workdir name ''''];
    [status,result] = system(cmd);
    disp(['Directory ' workdir name ' created']);
end



cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/inputs.tar.gz'''];
[status,result] = system(cmd);
if ( findstr(result,'No such file or directory'))
    disp(['Inputs are not present at brutus - creating locally ...']);          
%%------------------------------Creating the input file-------------------
[vars,output,perm_c] = CreateInputFiles_loop(opts,defopts);
AddtoSub('',0,0,0,0) %initialize
nr_batches = 0;
for i = 1:perm_c
   for j = 1:repeats
        batch = '';
        for k = 1:size(vars,1)
            if ( ~(strcmpi(char(output{k,i}),'Inf') || strcmpi(char(output{k,i}),'-Inf') )) %dont add any Inf Values - Fortran cant handel this via Input file
                                                                                            %have to handel this via the default values of pCMALib
                if (strcmpi(char(vars(k)),'output_folder') || strcmpi(char(vars(k)),'seed_folder'))
                    x = [char(vars(k)) ' = output' num2str(i) 'repeat' num2str(j)];
                else
                    if (strcmpi(char(vars(k)),'USE_LJC') )
                        x = ['USE_LJ = ' char(output{k,i})];
                    else
                        if (strcmpi(char(vars(k)),'STOPTIMEUSE') )
                            x = ['STOPTIME = ' char(output{k,i})];
                        else
                            x = [char(vars(k)) ' = ' char(output{k,i})];
                        end
                    end
                end
                batch = AddLine(batch,x);
            end
        end
        file = strcat('output',num2str(i),'repeat',num2str(j),'.txt');
        fid = fopen(file,'wt');
        fprintf(fid,'%s',batch);
        fclose(fid);
        % create the according submit batches
        nr_batches = AddtoSub(pathtoCMA,queuetime*60*60,expt,strcat('output',num2str(i),'repeat',num2str(j)),nrproc);
   end
end
%flush the final batch
nr_batches = AddtoSub(pathtoCMA,0,0,'flush',0);
%files are created - tar.gz them
cmd = ['tar zcvf inputs.tar.gz *output*repeat*'];
[status,result] = system(cmd);
cmd = ['tar zcvf batches.tar.gz *batch.txt'];
[status,result] = system(cmd);
cmd = ['rm *.txt'];
[status,result] = system(cmd);
%copy to brutus
cmd = ['scp batches.tar.gz brutus.ethz.ch:/' workdir name];
[status,result] = system(cmd);
cmd = ['scp inputs.tar.gz brutus.ethz.ch:/' workdir name];
[status,result] = system(cmd);
%extract the batches
cmd = ['ssh brutus.ethz.ch ''tar -xf ' workdir name '/batches.tar.gz -C ' workdir name '/'''];
[status,result] = system(cmd);
end %end input files not present


%% now input file should be there - check if we should submit
cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/*batch.txt | wc -w'''];
[status,result] = system(cmd);

if (isempty(findstr('No such file', result)))
batchnumber = str2num(result);
    
cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/sub* | wc -w'''];
[status,result] = system(cmd);

if (~isempty(findstr('No such file', result))) %we didnt submit yet
    %create the submit script
    runme = AddLine('','#!/bin/sh');
    runme = AddLine(runme, ['for n in {0..' num2str(batchnumber-1) '}; do']);
    if (nrproc < 2)
        runme = AddLine(runme,['bsub -W ' num2str(queuetime) ':00 -o report${n} < ${n}batch.txt > sub${n}']);
    else
        runme = AddLine(runme,['bsub -n ' num2str(nrproc) ' -W ' num2str(queuetime) ':00 -o report${n} < ${n}batch.txt > sub${n}']);
    end
    runme = AddLine(runme,['done']);
    
    fid = fopen('runme','wt');
    fprintf(fid,'%s',runme);
    fclose(fid);
    cmd = ['scp runme brutus.ethz.ch:/' workdir name];
    [status,result] = system(cmd);
    cmd = ['ssh brutus.ethz.ch ''chmod 777 ' workdir name '/runme'''];
    [status,result] = system(cmd);
    disp('Submitting Jobs - this may take quite a while');
    cmd = ['ssh brutus.ethz.ch ''source /etc/profile; cd ' workdir name '/; ./runme'''];
    [status,result] = system(cmd);
    disp('Finished submitting - call the script again later to gather the results');
    return;
else %there are already job submitted 
    disp('Jobs seem to be already submitted');
end

%jobs are already submitted
subnumber = str2num(result);
if (subnumber ~= batchnumber)
   disp('Something is wrong - the number of batches and submit differs - please check manualy');  
   return;
end

%%check if the results are in yet
cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/report* | wc -w'''];
[status,result] = system(cmd);
if (~isempty(findstr('No such file', result)))
    disp('Jobs seem submitted but none is yet finished - please check again later');
    return
else
    reportnumber = str2num(result);
    if (reportnumber ~= subnumber)
        disp([num2str(reportnumber) ' of ' num2str(subnumber) ' jobs are finished - please check again later']);
        return
    else
        cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/db_done | wc -w'''];
        [status,result] = system(cmd);
        if (~isempty(findstr('No such file', result)))
            disp('all jobs are finished - moving results to DB');
            cmd = ['ssh brutus.ethz.ch ''cd ' workdir name '/; tar zcvf my_outputs.tar.gz output*'''];
            [status,result] = system(cmd);
            cmd = ['scp brutus.ethz.ch:/' workdir name '/my_outputs.tar.gz ./'];
            [status,result] = system(cmd);
            cmd = ['tar -xf my_outputs.tar.gz'];
            [status,result] = system(cmd);
            cmd = ['ssh brutus.ethz.ch ''ls ' workdir name '/output*.tar.gz | wc -w'''];
            [status,result] = system(cmd);
            totalnumber = str2num(result);        
            FillDB(name,totalnumber/repeats, repeats,dbname,nrproc);
            cmd = ['ssh brutus.ethz.ch ''echo db_done >' workdir name '/db_done '''];
            [status,result] = system(cmd);
        else
            disp('Jobs are already moved to DB - delete db_done in the work dir to copy anyway');
        end
        
    end
end
else
    disp('Something is wrong - there are no batch files on brutus!');
end
end % end function Brutus_Benchmarks_final

