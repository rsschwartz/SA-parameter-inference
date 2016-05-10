function [ y ] = assemblyRMSD ( jobList, dataList, paramFile )

    load(paramFile);
    jobNum = length(jobList);
    curveList = cell(1,jobNum);
    jobID = cell(1,jobNum);
    %order = ceil(log(jobNum+1)/log(10));
    %formatStr = [repmat('\b',1,order),'%',num2str(order),'d'];
    formatStr = ['%',num2str(ceil(log(jobNum+1)/log(10))),'d'];
    submitJob = zeros(1,jobNum);
    activeJob = zeros(1,jobNum);
    finishJob = false(1,jobNum);
    indexList = 1 : jobNum;
    maximer = inf;
    %activeNum = 0;
    
    if ~isfield(qInfo,'interval'), qInfo.interval = 10; end
    
    if ~isfield(qInfo,'maxID'), qInfo.maxID = 10; end
    
    if ~isfield(qInfo,'maxJob'), qInfo.maxJob = jobNum; end
    
    disp('Submitting jobs...');

    while any(~finishJob)
        
        tic;
        
        if isfield(fileInfo,'debugFile2') && ...
                exist(fileInfo.debugFile2,'file')
            
            keyboard;
            
        end
        
        while sum(activeJob~=0) >= sum(qInfo.maxJob)
            
            pause(qInfo.interval);
            activeJob = countJob(activeJob,jobID,qInfo,fileInfo);
            
        end
        
        if any(submitJob==0)
        
            for q = 1 : numel(qInfo.submit)
                
                if sum(activeJob==q) < qInfo.maxJob(q), break, end
                
            end
            
            i = find(submitJob==0,1);
            jobScript = jobList{i};
            [status,ID] = system([qInfo.submit{q},' ',jobScript]);
        
            if status ~= 0
            
                disp('Error submitting job!');
                disp([qInfo.submit,' ',jobScript]);
                
            else
                
                submitJob(i) = q;
                jobID{i} = ID(1:find(ID=='.',1)-1);
                %{
                if length(ID) > qInfo.maxID, jobID{i} = ID(1:qInfo.maxID);
            
                else jobID{i} = ID;
                    
                end
                %}
            
            end
        
        end
        
        
        activeJob = countJob(submitJob.*~finishJob,jobID,qInfo,fileInfo);
        fprintf(['\rTotal: %d; submitted: %d; running: ',formatStr],...
            jobNum,sum(submitJob~=0),sum(activeJob~=0));
        
        for j = indexList(submitJob~=0 & ~finishJob & activeJob==0)
            
            currCurve = importSim(dataList{j});
            nmer = size(currCurve,2) - 1;
            
            if ~isempty(currCurve) && ...
                    isreal(currCurve) && ...
                    ~isnan(sum(sum(currCurve))) && ...
                    (~isfield(simInfo,'maximer') || nmer >= simInfo.maximer) 
                
                curveList{j} = currCurve;
                %curveList{j}(:,2:end) = currCurve(:,2:end) * diag(1:nmer)...
                %    / sum(currCurve(1,2:end));
                finishJob(j) = true;
                maximer = min(maximer,nmer);
                
            else
                
                submitJob(j) = false;
                
            end
            
        end
        
        pause(qInfo.interval - toc);
        
    end
    
    fprintf('\nAll jobs finished.\nCalculating objectives...\n');
    concNum = numel(simInfo.conc);
    tmpScore = zeros(concNum,1);
    simNum = jobNum / simInfo.repeat / concNum;
    y = zeros(simNum,1);
    
    for s = 1 : simNum
            
        for i = 1 : concNum
                
            [~,asmCurve] = avgCurve([],simInfo.data{i}(:,1),(1:maximer),...
                [],[],[],curveList(1:simInfo.repeat));
            asmCurve = asmCurve * diag(1:maximer) / sum(asmCurve(1,:));
            asmCurve = asmCurve - simInfo.data{i}(:,2:end);
            tmpScore(i) = norm(asmCurve,'fro')^2 / size(simInfo.data{i},1);
            curveList(1:simInfo.repeat) = [];
                %should further divide tmpScore(i) by
                % size(simInfo.data{i},2)
                %to make it normalizable
        end
        
        y(s) = (mean(tmpScore))^0.5;
            
    end
    
end

function [ activeJob ] = countJob ( activeJob, jobID, qInfo, fileInfo )

    jobIndex = 1 : numel(activeJob);
    [status,qStatus] = system(qInfo.query);
    
    while status ~= 0 && exist('qStatus','var') && ~isempty(qStatus)
        
        disp('Error querying queue!');
        disp(qStatus);
        disp(['Retry in ',num2str(qInfo.interval),' seconds...']);
        pause(qInfo.interval);
        clear qStatus;
        [status,qStatus] = system(qInfo.query);
        
        if isfield(fileInfo,'debugFile3') && ...
                exist(fileInfo.debugFile3,'file')
            
            keyboard;
            
        end
        
    end
    
    for i = jobIndex(activeJob~=0)
        
        if isempty(strfind(qStatus,jobID{i})), activeJob(i) = 0; end;
        
    end
        
end
