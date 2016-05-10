function [ simCurve, asmCurve, curveList ] = avgCurve ( simList, ...
    timeList, asm, K, timeStep, curveList, asmList )

    if ~exist('timeStep','var'), timeStep = 1; end
    
    if numel(timeList) == 1, timeList = (0:timeStep:timeList)'; end
    
    if ~exist('asm','var'), asm = []; end
    
    if ~exist('K','var') || isempty(K), K = 1; end
    
    if exist('curveList','var') && ~isempty(curveList)
        
        simNum = numel(curveList);
        
    else
        
        if ~isempty(simList)
        
            simNum = numel(simList);
        
        elseif ~isempty(asmList)
        
            simNum = numel(asmList);
            
        end
        
        curveList = cell(simNum,1);
        
    end
    
    timePoint = length(timeList);
    simCurve = zeros(timePoint,simNum);
    asmCurve = zeros(timePoint,length(asm),simNum);
    
    for i = 1 : simNum

        if ~isempty(curveList{i})
            
            currCurve = curveList{i};
            simPoint = size(currCurve,1);
            asmData = zeros(simPoint+1,1);
            
        else
            
            if ~isempty(simList)
                
                simData = importSim(simList{i});
            
                if sum(simData(end,2:end-1)) == 0 && ...
                        simData(end,1) < timeList(end)
               
                    simData = [simData;[timeList(end),simData(end,2:end)]];
                
                end
                
            elseif ~isempty(asmList)
            
                simData = asmList{i};
                
            end
            
            simPoint = size(simData,1);
            currCurve = zeros(simPoint,2);
            currCurve(:,1) = simData(:,1);
            simData(:,1) = [];
            asmData = simData(:,asm);
            mer = (1:size(simData,2))';
            currCurve(:,2) = simData*(mer.^2)./(simData*mer);
            curveList{i} = currCurve;
            
        end
        
        if currCurve(1,1) ~= 0
            
            currCurve = [[0,1];currCurve];
            simPoint = simPoint + 1;
        
        end
        
        count = 1;
        
        for j = 1 : timePoint;
            
            simTime = currCurve(count,1);
            stdTime = timeList(j);
            
            while simTime < stdTime && count < simPoint
                
                count = count + 1;
                simTime = currCurve(count,1);
                
            end

            if simTime == stdTime || ...
                    currCurve(count,1) == currCurve(count-1,1)
            
                simCurve(j,i) = currCurve(count,2);
                asmCurve(j,:,i) = asmData(count,:);
                
            else
                
                simCurve(j,i) = (currCurve(count,2) - currCurve(count-1,2)) / ...
                    (currCurve(count,1) - currCurve(count-1,1)) * (stdTime - ...
                    currCurve(count-1,1)) + currCurve(count-1,2);
                asmCurve(j,:,i) = (asmData(count,:) - asmData(count-1,:)) / ...
                    (currCurve(count,1) - currCurve(count-1,1)) * (stdTime - ...
                    currCurve(count-1,1)) + asmData(count-1,:);
                
            end
            
        end
        
    end
    
    simCurve = K*mean(simCurve,2);
    simCurve = [timeList,simCurve];
    asmCurve = mean(asmCurve,3);
    
end
