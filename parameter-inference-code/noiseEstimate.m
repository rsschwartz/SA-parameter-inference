function [ sa, ka, ca ] = noiseEstimate ( simList, simInfo, repeat )

    concNum = length(simInfo.conc);
    simNum = length(simList) / concNum;
    sa = zeros(1,repeat);
    ka = zeros(1,repeat);
    ca = cell(1,repeat);
    
    for i = 1 : repeat
        
        fprintf('Iteration %d out of %d\n',i,repeat);
        currList = cell(simInfo.repeat*concNum,1);
        currIndex = randi(simNum,1,simInfo.repeat);
        %currIndex = randperm(simNum);
        %currIndex = currIndex(1:simInfo.repeat);
        %unique(randi(simNum,1,simInfo.repeat));
        
        %while numel(currIndex) < simInfo.repeat
        
        %    currIndex = unique([currIndex,randi(simNum,1,...
        %        simInfo.repeat-numel(currIndex))]);
        
        %end
    
        for j = 1 : concNum
        
            currList(simInfo.repeat*(j-1)+1:simInfo.repeat*j) = ...
                simList(simNum*(j-1)+currIndex);
            
        end
        
        [sa(i),ka(i),ca{i}] = fitMultiCurve(currList,simInfo);
        fprintf('Current score: %d\n\n',sa(i));

    end
       
end