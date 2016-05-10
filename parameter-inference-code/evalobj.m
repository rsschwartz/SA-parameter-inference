function [ rmsd ] = evalobj ( paramFile, offset )

% calculate RMSD w.r.t. given offset. Additional parameters are passed over
% by paramFile.

    load(paramFile);
    
    if ~isfield(optInfo,'BS'), optInfo.BS = []; end

    if ~isfield(optInfo,'CS'), optInfo.CS = []; end
    
    if isfield(fileInfo,'log') && ~isempty(fileInfo.log)
        
        logFile = fopen(fileInfo.log,'a');
        
    end
    
    if ~isfield(optInfo,'iter'), optInfo.iter = 0; end
    
    if ~isfield(optInfo,'offset'), optInfo.offset = []; end
    
    if ~isfield(optInfo,'score'), optInfo.score = []; end
    
    if ~isfield(fileInfo,'log') || isempty(fileInfo.log) || logFile < 0
        
        logFile = fopen('/dev/null');
    
    end
    
    if ~isfield(fileInfo,'prefix'), fileInfo.prefix = []; end
    
    if isfield(fileInfo,'debugFile1') && ...
            exist(fileInfo.debugFile1,'file')
    
        keyboard;
    
    end
    
    prefix = blanks(20*numel(offset));
    ptr = 1;
    
    for i = 1 : numel(offset);
        
        count = numel(num2str(offset(i)));
        prefix(ptr:ptr+count) = [num2str(offset(i)),'_'];
        ptr = ptr + count + 1;
        
    end
    
    prefix(ptr-1:end) = [];
    optInfo.iter = optInfo.iter + 1;
    prefix = [prefix,'_',num2str(optInfo.iter),'_'];
    
    [bs,cs] = setParameter(optInfo.BS,optInfo.CS,optInfo.paraList,offset);
    jobPrefix = [fileInfo.prefix,'_',prefix];
    [jobList,dataList] = prepareJob(jobPrefix,[],[],bs,cs,fileInfo,simInfo);
    %sendJob(jobList,qInfo);
    %rmsd = fitMultiCurve(dataList,simInfo);
    
    rmsd = optInfo.funobj(jobList,dataList,paramFile);
    optInfo.offset = [optInfo.offset;offset];
    optInfo.score = [optInfo.score;rmsd];
    save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');
    
    fprintf(logFile,'Iteration: %d\nOffset: %s\nScore: %e\n\n',...
        optInfo.iter,prefix,rmsd);
    fclose(logFile);
    cleanup(fileInfo,jobPrefix);
    
end
