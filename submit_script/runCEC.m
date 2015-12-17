cmapath = '/cluster/home/infk/ofgeorg/checkout9/bin/libpcma';
workdir = '/cluster/work/infk/ofgeorg/';
repeats = 25;
queuetime = 1;
%Brutus_Benchmarks_final(cmapath,workdir,'testingname4',repeats,myopts,queuetime,'merge2');
%Brutus_Benchmarks_final(cmapath,workdir,'testingname5',repeats,myopts,queuetime,'merge2');
%Brutus_Benchmarks_final(cmapath,workdir,'ljtest',repeats,myopts,1,'ljtest'
%)
dimensions = [10,30,50];
uBound = [100.0, 100.0, 100.0, 100.0, 100.0, ...
                   100.0, 600.0, 32.0, 5.0, 5.0,...
                   0.5,...
                   3.14159265358979323846264338327950288419716939937510,...
                   1.0, 100.0, 5.0,...
                   5.0, 5.0, 5.0, 5.0, 5.0,...
                   5.0, 5.0, 5.0, 5.0, 5.0];
lBound = [-100.0, -100.0, -100.0, -100.0, -100.0, ...
                   -100.0, 0.0, -32.0, -5.0, -5.0,...
                   -0.5, ...
                   -3.14159265358979323846264338327950288419716939937510,...
                    -3.0, -100.0, -5.0,...
                   -5.0, -5.0, -5.0, -5.0, -5.0,...
                   -5.0, -5.0, -5.0, -5.0, 2.0 ];
               
g_min = [-4.5000000e+002 -4.5000000e+002 -4.5000000e+002 -4.5000000e+002 -3.1000000e+002 ...
          3.9000000e+002 -1.8000000e+002 -1.4000000e+002 -3.3000000e+002 -3.3000000e+002 ...
          9.0000000e+001 -4.6000000e+002 -1.3000000e+002 -3.0000000e+002  1.2000000e+002 ...
          1.2000000e+002  1.2000000e+002  1.0000000e+001  1.0000000e+001  1.0000000e+001 ...
          3.6000000e+002  3.6000000e+002  3.6000000e+002  2.6000000e+002  2.6000000e+002];
acc_rec = [1e-6 1e-6 1e-6 1e-6 1e-6 ...
           1e-2 1e-2 1e-2 1e-2 1e-2 ...
           1e-2 1e-2 1e-2 1e-2 1e-2 ...
           1e-2 1e-1 1e-1 1e-1 1e-1 ...
           1e-1 1e-1 1e-1 1e-1 1e-1];
dimensions = [10,30,50];
       
for (dimi = 1:size(dimensions,2))       
    for fcnr = 1:25       
    myopts = Brutus_Benchmarks_final();       
    myopts.benchmark = 'true';
    myopts.use_CEC = 'true';
    myopts.Benchfctnr = num2str(fcnr);
    myopts.CECFolders = '/cluster/home/infk/ofgeorg/checkout9/bin/';
    myopts.dimensions = num2str(dimensions(dimi));
    myopts.global_min = num2str(g_min(fcnr));
    myopts.accuracy = '1.E-8';
    myopts.record_accuracy = num2str(acc_rec(fcnr));
    myopts.restart_cma = 'true';
    myopts.restart_type = '1';
    myopts.MaxIncFac = '20';
    myopts.IncPopSize = '1.25';
    myopts.rel_sigma = '0.2';
    if (fcnr == 7 || fcnr == 25) 
        myopts.use_init_bounds = 'true';
        myopts.init_lBounds = num2str(lBound(fcnr));
        myopts.init_uBounds = num2str(uBound(fcnr));    
    else
        myopts.alldim_lBounds = num2str(lBound(fcnr));
        myopts.alldim_uBounds = num2str(uBound(fcnr));    
    end
    myopts.StopMaxFunEvals = num2str(dimensions(dimi) * 10000);
    if (dimensions(dimi) == 30)
        myopts.record_besthist = 'true';
        myopts.record_modulo = '1000';
    end
    expt = CEC2005runtime(dimi,fcnr);
    Brutus_Benchmarks_final(cmapath,workdir,['CEC_dim' num2str(dimi) 'fcnr' num2str(fcnr)],repeats,myopts,queuetime,'CEC',expt);
    end
end

