function [ scoreList, K, simCurve ] = fitMultiCurve ( simList, simInfo,...
    curveList )

    concNum = numel(simInfo.conc);
    scoreTemp = zeros(concNum,1);
    expCurve = simInfo.data;
    expPoint = zeros(concNum,1);
    expTime = cell(concNum,1);
    %simCurve = cell(concNum,1);
    
    for i = 1 : concNum
%{        
        if ~knownK
        
            expCurve{i}(:,2) = expCurve{i}(:,2) - ...
                mean(expCurve{i}(1:simInfo.lagPoint(i),2));
            
        end
%}        
        expPoint(i) = size(expCurve{i},1);
        expTime{i} = expCurve{i}(:,1);
        %simCurve{i} = zeros(expPoint(i),simNum);
        
    end
    
    if exist('curveList','var') && ~isempty(curveList)
        
        simNum = numel(curveList) / simInfo.repeat / concNum;
        simCurve = cell(simNum,concNum);
        
        for s = 1 : simNum
            
            for i = 1 : concNum
                
                tempCurve = avgCurve([],expTime{i},[],[],[],...
                    curveList(1:simInfo.repeat));
                simCurve{s,i} = tempCurve(:,2);
                curveList(1:simInfo.repeat) = [];
                
            end
            
        end
        
    elseif ~isempty(simList)
        
        simNum = numel(simList) / simInfo.repeat / concNum;
        simCurve = cell(simNum,concNum);
        
        for s = 1 : simNum
        
            for i = 1 : concNum
            
                tempCurve = avgCurve(simList(1:simInfo.repeat),expTime{i});
                simCurve{s,i} = tempCurve(:,2);
                simList(1:simInfo.repeat) = [];
                    
            end
                
        end
            
    else
            
        disp('Error: cannot find simulated data.');
        return
        
    end
       
    %simNum = numel(simList)/simInfo.repeat/concNum;
    scoreList = zeros(simNum,1);
    %curList = simList;
    
    if ~isfield(simInfo,'lagPoint') || isempty(simInfo.lagPoint)
        
        simInfo.lagPoint = ones(concNum,1);
        
    elseif numel(simInfo.lagPoint) == 1 && concNum > 1
        
        simInfo.lagPoint = ones(concNum,1) * simInfo.lagPoint;
        
    end
    
    if isfield(simInfo,'K') && ~isempty(simInfo.K)
            
        K = ones(simNum,1) * simInfo.K;
        
    else
        
        K = -ones(simNum,1);
            
    end
    
    if K(1) > 0, knownK = true;
        
    else knownK = false;
        
    end
    
    if ~knownK
        
        for i = 1 : concNum
            
            %if expCurve{i}(1,1) == 0
        
                expCurve{i}(:,2) = expCurve{i}(:,2) - ...
                    mean(expCurve{i}(1:simInfo.lagPoint(i),2));
                
            %end
            
        end
        
        %expPoint(i) = size(expCurve{i},1);
        %expTime{i} = expCurve{i}(:,1);
        %simCurve{i} = zeros(expPoint(i),simNum);
        
    end
    
    for s = 1 : simNum
%{        
        for i = 1 : concNum
            
            tempCurve = avgCurve(curList(1:simInfo.repeat),expTime{i});
            simCurve{i}(:,s) = tempCurve(:,2);
            curList(1:simInfo.repeat)=[];
            
        end
%}        
        if size(K,2) == 1 &&  K(s) <= 0
            
            Knumer = zeros(concNum,1);
            Kdenom = zeros(concNum,1);
            
            for i = 1 : concNum
                
                %Knumer(i) = expCurve{i}(:,2)'*(simCurve{s,i}-simCurve{s,i}(1)) * ...
                %    simInfo.conc(i) / expPoint(i);
                %Kdenom(i) = (simCurve{s,i}-simCurve{s,i}(1))'*(simCurve{s,i}-simCurve{s,i}(1)) *...
                %    simInfo.conc(i)^2 / expPoint(i);
                Knumer(i) = expCurve{i}(:,2)'*(simCurve{s,i}-1) * ...
                    simInfo.conc(i) / expPoint(i);
                Kdenom(i) = (simCurve{s,i}-1)'*(simCurve{s,i}-1) *...
                    simInfo.conc(i)^2 / expPoint(i);
            
            end
            
            if sum(Knumer+Kdenom) == 0
                
                K(s) = 1;
                
            else
                
                K(s) = sum(Knumer)/sum(Kdenom);
                
            end
            
        elseif size(K,2) == concNum && K(s,i) <= 0
            
            for i = 1 : concNum
            
                if sum(simCurve{i}(:,s)-1) == 0
                    
                    K(s,i) = 1;
                    
                else
                    
                    %K(s,i) = (expCurve{i}(:,2)'*(simCurve{s,i}-simCurve{s,i}(1))) / ...
                    %    ((simCurve{s,i}-simCurve{s,i}(1))'*(simCurve{s,i}-simCurve{s,i}(1))) / ...
                    %    simInfo.conc(i);
                    K(s,i) = (expCurve{i}(:,2)'*(simCurve{s,i}-1)) / ...
                        ((simCurve{s,i}-1)'*(simCurve{s,i}-1)) / ...
                        simInfo.conc(i);
                    
                end
                    
            end
            
        end
        
        for i = 1 : concNum
            
            if size(K,2) == 1
            
                if knownK
                    
                    diff = expCurve{i}(:,2) - ...
                        simCurve{s,i}*K(s)*simInfo.conc(i);
                    
                else
                    
                    diff = expCurve{i}(:,2)+K(s)*simInfo.conc(i)*simCurve{s,i}(1) - ...
                        simCurve{s,i}*K(s)*simInfo.conc(i);
                    
                end
                
            elseif size(K,2) == concNum
                
                if knownK
                
                    diff = expCurve{i}(:,2) - ...
                        simCurve{s,i}*K(s,i)*simInfo.conc(i);
                    
                else
                    
                    diff = expCurve{i}(:,2)+K(s,i)*simInfo.conc(i)*simCurve{s,i}(1) - ...
                        simCurve{s,i}*K(s,i)*simInfo.conc(i);
                    
                end
                
            end
            
            scoreTemp(i) = diff'*diff/expPoint(i);
            
        end
        
        if size(K,2) == 1

            scoreList(s) = (sum(scoreTemp)/concNum)^0.5;
        
        elseif size(K,2) == concNum

            scoreList(s) = sum(scoreTemp.^0.5)/concNum;

        end

    end
    
end
