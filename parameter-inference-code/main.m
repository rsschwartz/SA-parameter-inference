%% Input arguments
%optional items can be left blank or commented.

fileInfo.xmlRule = '/home/lxie1/simulator/hbv54.xml';
%fileInfo.xmlRule{2} = '/home/lxie1/simulator/hpv198.xml';
%fileInfo.xmlRule{3} = '/home/lxie1/simulator/hpv220.xml';
%where your xml file is

fileInfo.command = ['java -Xmx1024m -jar '...
    '/home/lxie1/simulator/dessa1.5.8.jar'];
%command of running simulator. You can replace it with C++ simulator.

fileInfo.workFolder = '/home/lxie1/hbvSLS2/';
%folder for storing job scripts, input and output files.

fileInfo.errFolder = [fileInfo.workFolder,'err/'];
%folder for error information

%fileInfo.preset = ['#PBS -o ',fileInfo.errFolder,...
%    ' -e ',fileInfo.errFolder,' -l pmem=1gb'];
%fileInfo.preset = '#PBS -o /dev/null -e /dev/null';
%preset running environment. eg. PBS stuff


fileInfo.tempFolder = '/scratch/';
%folder for storing temporary data files. OPTIONAL

fileInfo.log = [fileInfo.workFolder,'hbv3_140920.log'];
%log file of simulation progress.
%OPTIONAL, default = print to screen only; get your pen and paper ready

fileInfo.debugFile1 = [fileInfo.workFolder,'debug1'];
fileInfo.debugFile2 = [fileInfo.workFolder,'debug2'];
fileInfo.debugFile3 = [fileInfo.workFolder,'debug3'];
fileInfo.debugFile4 = [fileInfo.workFolder,'debug4'];
%debug files for different level of debugging.

qInfo.submit{1} = ['ssh -x lanec1 /opt/torque/bin/qsub',...
    ' -q rs1',...
    ' -o ',fileInfo.errFolder,...
    ' -e ',fileInfo.errFolder,...
    ' -l pmem=1gb,pvmem=2gb,walltime=02:00:00'];
qInfo.submit{2} = ['ssh -x lanec1 /opt/torque/bin/qsub',...
    ' -q pool1',...
    ' -o ',fileInfo.errFolder,...
    ' -e ',fileInfo.errFolder,...
    ' -l pmem=1gb,pvmem=2gb,walltime=02:00:00'];
%command for submitting jobs, including parameters.

qInfo.query = 'qstat|grep lxie1';
%command for querying queue status.

qInfo.interval = 0.1;
%time interval for checking queue status.
%OPTIONAL, default = 10;

qInfo.maxID = 10;
%the ID returned by qsub is really long; decide how many characters you
%would like to store.
%OPTIONAL, default = 10;

qInfo.maxJob = [90,60];
%maximum number of jobs hanging on the queue.
%OPTIONAL. Set it when your colleagues complain.

noise = 0.025; % real 3-SLS, 50 repeat
%background noise level for score calculation.

simInfo.time = 550;
%simulation time as an argument to DESSA.

simInfo.repeat = 50;
%how many times you want to repeat on one set of parameter.

%simInfo.predictRepeat = 250;

simInfo.interval = 1e4;
%events interval for DESSA output.

simInfo.maximer = 120;

%load /home/lxie1/simulator/sim_hbv3_f19b39_f25b45;
%simInfo.data = {shpv145,shpv198,shpv220};
load /home/lxie1/simulator/hbv;
simInfo.data = {hbv54, hbv82, hbv108};
%simInfo.data = {sim_hbv_sls_f19b39_f25b45};
%simInfo.data = {sls1_5,sls1,sls2};
%load /home/lxie1/simulator/hpv145;
%load /home/lxie1/simulator/hpv198;
%load /home/lxie1/simulator/hpv220;
%simInfo.data{1} = hpv145;
%simInfo.data{2} = hpv198;
%simInfo.data{3} = hpv220;
%experimental data extracted from the light scattering curve

simInfo.conc = [5.4; 8.2; 10.8];
%concentrations of coat proteins in experiment. The absolute unit was not
%important; but the unit must be consistent. Must be column vector.

simInfo.lagPoint = [100,50,20];
%the #points in the lag phase of each curve. Don't need to be accurate but
%fewer is better. Optional. Default = 1.

simInfo.K = [];
%scaling factor of Rayleigh ratio combined with concentration and molecular
%weight
%{
optInfo.BS = {'bst0a' 'bst0b' [1.963204407833257,-4.882236855678068]
              'bst0c' 'bst1a' [1.963204407833257,-4.952929196273922]
              'bst0d' 'bst1b' [1.963204407833257,-4.948210637772643]
              'bst1c' 'bst1d' [1.963204407833257,-5.519189689266716]};
%}
optInfo.BS = {'bst0a' 'bst0b' [log(91.8764758528154)/log(10),log(7.35886889953384e-06)/log(10)]
              'bst0c' 'bst1a' [log(91.8764758528154)/log(10),log(1.171150717264e-05)/log(10)]
              'bst0d' 'bst1b' [log(91.8764758528154)/log(10),log(1.171150717264e-05)/log(10)]
              'bst1c' 'bst1d' [log(91.8764758528154)/log(10),log(9.58960364746733e-06)/log(10)]};
%
%binding sites and log-10 based binding / breaking / fastbinding time
%you can just write the sites and times involved in optimization; you don't
%have to write all information defined in the xml file.
%OPTIONAL. You don't have to write this term if you don't want to play with
%binding times.

%optInfo.CS = {'bs0' 'bs1' 4
%              'bs1' 'bs0' 4};
%conformation switch time for optimization.
%OPTIONAL. You don't have to write this term if you don't want to play with
%conformation switch.

%optInfo.initBias = 1;
%initial bias. Bias equals to s in Senthil's paper and lambda in Russell's
%book.

%optInfo.maxBias = 1024;
%maximum Bias.

%optInfo.getS = @(s)1/s;
%get grid size from bias. In Senthil's paper grid size = 1/s.

optInfo.minGridSize = 1/1024;
%optimization search will terminate if grid size reaches this threshold.

optInfo.initGridSize = 1;
%initial grid size for generating the first parameter candidates.

optInfo.eps = 1e-3;
%hreshold for improvement.
%OPTIONAL, default = 1e-3;

optInfo.noise = noise / simInfo.repeat^0.5;
%noise level adjusted by the number of repeats;

optInfo.initScore = 2.263900e-01;
%set the score if you know it!

optInfo.maxIter = 10;
%maximum # of consecutive iteration without improvements;

optInfo.lb = -1;%-1
optInfo.ub = 1;%1
%bounds for times. The lb here bounds the lower limit of binding time to 1,
%the lower limit of break time to 1e-6. The ub here bounds the upper limit
%of binding time to 1e6, the upper limit of break time to 1.

%optInfo.funobj = @assemblyRMSD;
optInfo.funobj = @computeObjective;
%function to compute objective

optInfo.iter = 0;
%variable that stores iteration number

RandStream.setDefaultStream(RandStream('mt19937ar','seed',647456));
%set the random seed; OPTIONAL

%fileInfo.paramFile = [fileInfo.workFolder,'paramFile.mat'];

format long
%
%% First, 2-parameter search, all 4 bind times vs. all 4 break times.

fileInfo.prefix = 'hbv3_FR';

%indicate the working folder: /home/lxie1/hbvSim/
%every file created by this batch of optimization will start with hbv_FR*.*

optInfo.paraList{1} = {'b',(1:4),1};
optInfo.paraList{2} = {'b',(1:4),2};
%group of parameters. Each unit in paraList corresponds to the coordinates
%in optInfo.BS or CS. The first row of paraList corresponds to BS, and the 
%second row for CS. The columns correspond to parameter groups. Here only 
%BS is considered; if you want to add conformation switch time, write:
%paraList{2,1} = ...
save paramFile fileInfo qInfo simInfo optInfo;
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj','paramFile',lb,ub);
optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);

%[newTime,score] = LMopt ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = optimizer ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = snobopt ( fileInfo, qInfo, simInfo, optInfo );
%LMopt uses Levenburg-Marquardt methods; optimizer uses the method in
%Senthil's paper.

%optInfo.BS = newTime.BS;
%optInfo.CS = newTime.CS;

%cleanup(fileInfo,fileInfo.prefix);

%% Second, 2-parameter search, 1&4 break times vs. 2&3 break time.

fileInfo.prefix = 'hbv3_14_23';

optInfo.paraList{1} = {'b',[1,4],2};
optInfo.paraList{2} = {'b',[2,3],2};
optInfo.initScore = fbest;
%optInfo.getS = @(s)1/s;

save paramFile fileInfo qInfo simInfo optInfo;
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj','paramFile',lb,ub);
optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);

%[newTime,score] = LMopt ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = optimizer ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = snobopt ( fileInfo, qInfo, simInfo, optInfo );

%optInfo.BS = newTime.BS;
%optInfo.CS = newTime.CS;

%cleanup(fileInfo,fileInfo.prefix);
%
%% Third, 3-parameter search, 1 break time vs. 2&3 break time vs. 4 break time

fileInfo.prefix = 'hbv3_mcs_1_23_4';

optInfo.paraList{1} = {'b',1,2};
optInfo.paraList{2} = {'b',[2,3],2};
optInfo.paraList{3} = {'b',4,2};

optInfo.initScore = fbest;
fileInfo.paramFile = [fileInfo.workFolder,fileInfo.prefix,'_param.mat'];
save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj',fileInfo.paramFile,lb,ub);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);

%[newTime,score] = LMopt ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = optimizer ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = snobopt ( fileInfo, qInfo, simInfo, optInfo );

%optInfo.BS = newTime.BS;
%optInfo.CS = newTime.CS;

cleanup(fileInfo,fileInfo.prefix);
%
%% Fourth, 4-parameter search, each break time

fileInfo.prefix = 'hbv3_mcs_1234';

optInfo.paraList{1} = {'b',1,2};
optInfo.paraList{2} = {'b',2,2};
optInfo.paraList{3} = {'b',3,2};
optInfo.paraList{4} = {'b',4,2};

optInfo.initScore = fbest;
fileInfo.paramFile = [fileInfo.workFolder,fileInfo.prefix,'_param.mat'];
save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj',fileInfo.paramFile,lb,ub);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);

%newTime = LMopt ( fileInfo, qInfo, simInfo, optInfo );
%newTime = optimizer ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = snobopt ( fileInfo, qInfo, simInfo, optInfo );

%optInfo.BS = newTime.BS;
%optInfo.CS = newTime.CS;

cleanup(fileInfo,fileInfo.prefix);
%{
%% Fifth, 5-parameter search, each break time plus unique bind time

fileInfo.prefix = 'hbv3_all5';

optInfo.paraList{1} = {'b',1,2};
optInfo.paraList{2} = {'b',2,2};
optInfo.paraList{3} = {'b',3,2};
optInfo.paraList{4} = {'b',4,2};
optInfo.paraList{5} = {'b',(1:4),1};

optInfo.initScore = fbest;
%optInfo.getS = @(s)0.25/s;
save paramFile fileInfo qInfo simInfo optInfo;
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj','paramFile',lb,ub);
optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);

%newTime = LMopt ( fileInfo, qInfo, simInfo, optInfo );
%newTime = optimizer ( fileInfo, qInfo, simInfo, optInfo );
%[newTime,score] = snobopt ( fileInfo, qInfo, simInfo, optInfo );

%optInfo.BS = newTime.BS;
%optInfo.CS = newTime.CS;

%cleanup(fileInfo,fileInfo.prefix);
%
%% All parameter search
%
fileInfo.prefix = 'hbv_all8';

optInfo.paraList{1} = {'b',1,1};
optInfo.paraList{2} = {'b',2,1};
optInfo.paraList{3} = {'b',3,1};
optInfo.paraList{4} = {'b',4,1};
optInfo.paraList{5} = {'b',1,2};
optInfo.paraList{6} = {'b',2,2};
optInfo.paraList{7} = {'b',3,2};
optInfo.paraList{8} = {'b',4,2};

%optInfo.initScore = score;
%optInfo.getS = @(s)0.25/s;
save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');
lb = optInfo.lb * ones(numel(optInfo.paraList),1);
ub = optInfo.ub * ones(numel(optInfo.paraList),1);
[xbest,fbest] = mcs('evalobj',fileInfo.paramFile,lb,ub);

if fbest < optInfo.initScore
    
    optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);

else
    
    fbest = optInfo.initScore;
    
end

logFile = fopen(fileInfo.log,'a');
fprintf(logFile,'\n\n%s\n\nBest score: %e\n\n',fileInfo.prefix,fbest);
fclose(logFile);
optInfo.BS = setParameter(optInfo.BS,[],optInfo.paraList,xbest);
%}