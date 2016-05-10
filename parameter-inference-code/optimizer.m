function [ predictPara, minScore ] = optimizer ( paramFile )

%Use the method discribed in Senthil's paper, but here more grids are
%collected during each step. Here "r" equals to "s" in his paper, while "s"
%here refers to grid size.

    load(paramFile);
    %improvement = inf;
    r = optInfo.initBias;
    nPara = size(optInfo.paraList,2);
    quadFun = makeQuadFun(nPara);
    predictRepeat = 20*nPara^2;
 
    if ~isfield(optInfo,'eps'), optInfo.eps = 1e-3; end
    
    if ~isfield(optInfo,'BS'), optInfo.BS = []; end

    if ~isfield(optInfo,'CS'), optInfo.CS = []; end
    
    if ~isfield(optInfo,'offset'), optInfo.offset = []; end
    
    if ~isfield(optInfo,'score'), optInfo.score = []; end
    
    if isfield(fileInfo,'log') && ~isempty(fileInfo.log)
        
        logFile = fopen(fileInfo.log,'a');
        
    end
    
    if ~isfield(fileInfo,'log') || isempty(fileInfo.log) || logFile < 0
        
        logFile = fopen('/dev/null');
    
    end
%{    
    if isfield(fileInfo,'resultFile')
        
        param_log = true;
        
        if exist(fileInfo.resultFile,'file'), load(fileInfo.resultFile); end
        
        if ~exist('scoreList','var'), scoreList = []; end
        
        if ~exist('paramList','var'), paramList = []; end
        
    else
        
        param_log = false;
    
    end
%}    
    fprintf('Begin search process for %s...\n\n',fileInfo.prefix);
    fprintf(logFile,'Begin search process for %s...\n\n',fileInfo.prefix);

    if ~isfield(optInfo,'initScore') || isempty(optInfo.initScore)
        
        disp('Calculating initial score...');
        jobPrefix = [fileInfo.prefix,'_init'];
        simInfo.repeat = predictRepeat;
        save(paramFile,'fileInfo','qInfo','simInfo','optInfo');
        [jobList,dataList] = prepareJob(jobPrefix,[],[],optInfo.BS,...
            optInfo.CS,fileInfo,simInfo);
        %sendJob(jobList,qInfo);
        %optInfo.initScore = fitMultiCurve(dataList,simInfo);
        optInfo.initScore = optInfo.funobj(jobList,dataList,paramFile);
        cleanup(fileInfo,jobPrefix);
        disp('Initial score calculated.');
        
        %if param_log == true
            
            optInfo.score = optInfo.initScore;
            optInfo.offset = zeros(1,nPara);
            
        %end
        
    end
    
    fprintf('Initial score: %e\n\n',optInfo.initScore);
    fprintf(logFile,'Initial score: %e\n\n',optInfo.initScore);
    minScore = optInfo.initScore;
    count = 0;
    
    while r <= optInfo.maxBias % && improvement >= optInfo.eps
            
        if isfield(fileInfo,'debugFile1') && ...
                exist(fileInfo.debugFile1,'file')
            
            keyboard;
            
        end
        
        count = count + 1;
        s = optInfo.getS(r);
        x = prepareGrid(s,nPara);
        jobList = [];
        dataList = [];
        simInfo.repeat = simInfo.sampleRepeat;
        save(paramFile,'fileInfo','qInfo','simInfo','optInfo');
        iterPrefix = [fileInfo.prefix,'_iter',num2str(count)];
        fprintf('Search iteration: %d\nGrid size: %e\n',count,s);
        fprintf(logFile,'Search iteration: %d\nGrid size: %e\n',count,s);
        
        for i = 1 : size(x,1)

            if sum(abs(x(i,:))) == 0, continue; end

            fprintf('\rPreparing sample points: %d out of %d',i,size(x,1));
            [bs,cs] = setParameter(optInfo.BS,optInfo.CS,...
                optInfo.paraList,x(i,:));
            jobPrefix = [iterPrefix,'/grid',num2str(i)];
            [jobList,dataList] = prepareJob(jobPrefix,jobList,...
                dataList,bs,cs,fileInfo,simInfo);

        end

        %sendJob(jobList,qInfo);
        %y = fitMultiCurve(dataList,simInfo);
        y = optInfo.funobj(jobList,dataList,paramFile);
        disp('Making prediction...');
        offset = getOffset(x,[minScore;y],nPara,r,s,optInfo.lb,...
            optInfo.ub,quadFun);
        %fprintf(logFile,'Sample scores:\n');
        %fprintf(logFile,'%e\t',[minScore;y]);
        fprintf(logFile,'Parameter offset:\n');
        fprintf(logFile,'%e\t',offset);
        cleanup(fileInfo,iterPrefix);
        
        simInfo.repeat = predictRepeat;
        save(paramFile,'fileInfo','qInfo','simInfo','optInfo');
        [bs,cs] = setParameter(optInfo.BS,optInfo.CS,optInfo.paraList,...
            offset);
        jobPrefix = [fileInfo.prefix,'_predict',num2str(count)];
        [jobList,dataList] = prepareJob(jobPrefix,[],[],bs,cs,...
            fileInfo,simInfo);
        %sendJob(jobList,qInfo);
        %newScore = fitMultiCurve(dataList,simInfo);
        newScore = optInfo.funobj(jobList,dataList,paramFile);
        fprintf('Score of new prediction: %e\n',newScore);
        fprintf(logFile,'\nScore of new prediction: %e\n',newScore);
        %cleanup(fileInfo,jobPrefix);
        
        %if param_log == true
            
            optInfo.score = [optInfo.score;y;newScore];
            optInfo.offset = [optInfo.offset;x;offset'];
            save(paramFile,'fileInfo','qInfo','simInfo','optInfo');
            
        %end
        
        if newScore < minScore

            %improvement = (minScore - newScore) / minScore;
            minScore = newScore;
            
            if ~isfield(optInfo,'minBias') || r >= 2*optInfo.minBias
                
                r = r/2;
            
            end
            
            optInfo.BS = bs;
            optInfo.CS = cs;
            disp('Improved!');
            fprintf(logFile,'Improved! New parameters:\n');
            
            for i = 1 : size(bs,1)
                
                fprintf(logFile,'%s %s %.15e %.15e\n',bs{i,1},bs{i,2},...
                    bs{i,3}(1),bs{i,3}(2));
                
            end
            
            for i = 1 : size(cs,1)
                
                fprintf(logFile,'%s %s %.15e\n',cs{i,1},cs{1,2},cs{i,3});
                
            end
            
        else

            disp('No improvement. Increase bias...');
            %improvement = inf;
            r = 2*r;

        end

        fprintf('Iteration %d finished.\n\n',count);
        fprintf(logFile,'Iteration %d finished.\n\n',count);
        
    end

    %if improvement < optInfo.eps
        
    %    disp('Minimal improvement reached.');
        
    %end
    
    if r > optInfo.maxBias
        
        disp('Minimal grid size reached.');
        
    end

    fprintf('End search process for %s...\n\n',fileInfo.prefix);
    fprintf(logFile,'End search process for %s...\n\n',fileInfo.prefix);
    fclose(logFile);
    predictPara.BS = optInfo.BS;
    predictPara.CS = optInfo.CS;
    
end

function [ x ] = prepareGrid( s, nPara )

    if nPara == 1, x = [0;s;-s]; return, end
    
    %x = zeros(1+2*nPara+4*nchoosek(nPara,2),nPara);
    %x = zeros(1+nPara+nPara^2,nPara);
    x = zeros(1+2*nPara^2,nPara);
    
    for i = 1 : nPara
        
        x(i*2,i) = s;
        x(i*2+1,i) = -s;
        
    end
    
    count = 1+2*nPara;
    
    for i = 1 : nPara
        
        for j = i + 1 : nPara
            
            x(count+1,[i,j]) = [s,s];
            x(count+2,[i,j]) = [-s,-s];
            x(count+3,[i,j]) = [s,-s];
            x(count+4,[i,j]) = [-s,s];
            count = count + 4;
            %count = count + 2;
            
        end
        
    end
       
end

function [ offset ] = getOffset (x,y,nPara,r,s,lb,ub,quadFun)

    c = nlinfit(x,y,quadFun,ones(1,(nPara+1)*(nPara+2)/2));
    lb = lb*s*ones(1,nPara);
    ub = ub*s*ones(1,nPara);
    quadFunFix = makeQuadFun(nPara,c);
    options = optimset('TolFun',min(y)*1e-6,'TolX',s*1e-6,'TolCon',s*1e-6);
    %rsOffset = simulannealbnd(quadFunFix,zeros(1,nPara),lb,ub)';
    rsOffset = fmincon(quadFunFix,zeros(1,nPara),[],[],[],[],lb,ub,[],...
        options)';
    G = (y(2:2:nPara*2) - y(3:2:1+nPara*2)) / (2*s);
    %H = getHessian(y,s,nPara);
    
    if norm(G) ~= 0, G = G / norm(G); end
        
    offset = (rsOffset-r*s*G) / (1+r);
    %offset = -(H + r*diag(diag(H))) \ G;
    %offset = (rsOffset-r*H\G) / (1+r);

end
%{
function [ H ] = getHessian ( y, s, nPara )

    if nPara == 1, H = (y(2)+y(3)-2*y(1))/(s^2); return, end;
    
    H = zeros(nPara,nPara);
    
    for i = 1 : nPara, H(i,i) = (y(i*4)+y(i*4+1)-2*y(1))/(s^2); end
    
    count = 1 + 2*nPara;
    
    for i = 1 : nPara
        
        for j = i+1 : nPara
            
            H(i,j) = (y(count+1)+y(count+2)-2*y(1))/(2*s^2) - ...
                (H(i,i)+H(j,j))/2;
            H(j,i) = H(i,j);
            count = count + 4;
            
        end
        
    end
    
end
%}
function [ quadFun ] = makeQuadFun ( nVar, coef )

    nCoef = (nVar+2)*(nVar+1)/2;
    coefStr = cell(1,nCoef);
    funStr = cell(1,nCoef+1);

    if nargin == 2 && length(coef) >= nCoef && isnumeric(coef)
    
        funStr(1) = {'@(x)'};
        
        for i = 1 : nCoef
            
            if coef(i) >= 0, coefStr(i) = {['+',num2str(coef(i))]};
            
            else coefStr(i) = {num2str(coef(i))};
                
            end
                
        end
        
    else

        funStr(1) = {'@(c,x)'};
        
        for i = 1 : nCoef, coefStr(i) = {['+c(',num2str(i),')']}; end
        
    end

    varStr = cell(1,nVar+1);
    varStr(1) = {'1'};
    
    for i = 1 : nVar, varStr(i+1) = {['x(:,',num2str(i),')']}; end
    
    term = 0;
    
    for i = 1 : nVar + 1
        
        for j = i : nVar + 1
            
            term = term + 1;
            c = coefStr{term};
            v1 = varStr{i};
            v2 = varStr{j};
            funStr(term+1) = {[c,'*',v1,'.*',v2]};
            
        end
        
    end
    
    quadFun = eval(cell2mat(funStr));

end
