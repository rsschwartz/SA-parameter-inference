function [ y ] = computeObjective ( jobList, dataList, paramFile )

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
    %activeNum = 0;
    
    if ~isfield(qInfo,'interval'), qInfo.interval = 10; end
    
    if ~isfield(qInfo,'maxID'), qInfo.maxID = 10; end
    
    if ~isfield(qInfo,'maxJob'), qInfo.maxJob = jobNum; end
    
    disp('Submitting jobs...');
%{    
    for i = 1 : min(jobNum,qInfo.maxJob)
       
        jobScript = jobList{i};
        [status,ID] = system([qInfo.submit,' ',jobScript]);
        
        while status ~= 0
            
            disp('Error submitting job!');
            disp([qInfo.submit,' ',jobScript]);
            disp(['Retry in ',num2str(qInfo.interval),' seconds...']);
            pause(qInfo.interval);
            [status,ID] = system([qInfo.submit,' ',jobScript]);
            
        end

        if length(ID) > qInfo.maxID, jobID{i} = ID(1:qInfo.maxID);
            
        else jobID{i} = ID;
            
        end
        
        activeJob(i) = true;
        fprintf('\rTotal: %d; submitted: %d; running: %s',jobNum,i,...
            blanks(order));
        pause(qInfo.interval);
        
    end
    
    activeJob = countJob(activeJob,jobID,qInfo);
    activeNum = sum(activeJob);
    fprintf(formatStr,activeNum);
%}        
    while any(~finishJob)
        
        tic;
        
        if isfield(fileInfo,'debugFile2') && ...
                exist(fileInfo.debugFile2,'file')
            
            keyboard;
            
        end
        
        while sum(activeJob~=0) >= sum(qInfo.maxJob)
            
            pause(qInfo.interval);
            activeJob = countJob(activeJob,jobID,qInfo,fileInfo);
            %activeNum = sum(activeJob);
            %fprintf(formatStr,activeNum);
            
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
                disp([qInfo.submit{q},' ',jobScript]);
                %disp(['Retry in ',num2str(qInfo.interval),' seconds...']);
                %pause(qInfo.interval);
                %[status,ID] = system([qInfo.submit,' ',jobScript]);
                
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
        
        %activeJob(i) = true;
        %activeNum = activeNum + 1;
        %fprintf('\rTotal: %d; submitted: %d; running: %s',jobNum,i,...
        %    blanks(order));
        %fprintf(formatStr,activeNum);
        
        activeJob = countJob(submitJob.*~finishJob,jobID,qInfo,fileInfo);
        %activeNum = sum(activeJob);
        fprintf(['\rTotal: %d; submitted: %d; running: ',formatStr],...
            jobNum,sum(submitJob~=0),sum(activeJob~=0));
        
        for j = indexList(submitJob~=0 & ~finishJob & activeJob==0)
            
            currCurve = importSim(dataList{j});
            
            if ~isempty(currCurve) && ...
                    isreal(currCurve) && ...
                    ~isnan(sum(sum(currCurve))) && ...
                    (~isfield(simInfo,'maximer') || ...
                    size(currCurve,2) >= simInfo.maximer + 1) 
                
                curveList{j} = zeros(size(currCurve,1),2);
                curveList{j}(:,1) = currCurve(:,1);
                mer = (1 : size(currCurve,2)-1)';
                curveList{j}(:,2) = currCurve(:,2:end) * mer.^2 ./ ...
                    (currCurve(:,2:end) * mer);
                finishJob(j) = true;
                
            else
                
                submitJob(j) = false;
                
            end
            
        end
        
        pause(qInfo.interval - toc);
        
    end
%{    
    while activeNum >= 1
        
        pause(qInfo.interval);
        activeJob = countJob(activeJob,jobID,qInfo);
        activeNum = sum(activeJob);
        fprintf(formatStr,activeNum);
        
    end
%}    
    fprintf('\nAll jobs finished.\nCalculating objectives...\n');
    y = fitMultiCurve([],simInfo,curveList);
    
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
