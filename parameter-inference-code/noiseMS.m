function [ score ] = noiseMS ( simList, simInfo, repeat )

    concNum = length(simInfo.conc);
    simNum = length(simList) / concNum;
    score = zeros(1,repeat);
    tmpScore = zeros(concNum,1);
    maximer = size(simInfo.data{1},2) - 1;
    
    for i = 1 : repeat
        
        fprintf('Iteration %d out of %d\n',i,repeat);
        currIndex = randi(simNum,1,simInfo.repeat);
        
        %{        
        currIndex = unique(randi(simNum,1,simInfo.repeat));
        
        while numel(currIndex) < simInfo.repeat
        
            currIndex = unique([currIndex,randi(simNum,1,...
                simInfo.repeat-numel(currIndex))]);
        
        end
        %}
        
        for j = 1 : concNum
        
            [~,asmCurve] = avgCurve(simList(simNum*(j-1)+currIndex),...
                simInfo.data{j}(:,1),(1:maximer));
            asmCurve = asmCurve * diag(1:maximer) / sum(asmCurve(1,:));
            asmCurve = asmCurve - simInfo.data{j}(:,2:end);
            tmpScore(j) = norm(asmCurve,'fro')^2 / size(simInfo.data{j},1);
                
        end
        
        score(i) = (mean(tmpScore))^0.5;
        fprintf('Current score: %d\n\n',score(i));

    end
       
end