function [ predictPara, minScore ] = snobopt ( paramFile )

    load(paramFile);
    tempFile = [fileInfo.workFolder,fileInfo.prefix,'.mat'];
    nPara = size(optInfo.paraList,2);
    
    if exist(tempFile,'file'), system(['rm -f ',tempFile]); end
    
    if ~isfield(optInfo,'eps'), optInfo.eps = 1e-3; end
    
    if ~isfield(optInfo,'BS'), optInfo.BS = []; end

    if ~isfield(optInfo,'CS'), optInfo.CS = []; end
    
    if ~isfield(optInfo,'initScore'), optInfo.initScore = inf; end
    
    if ~isfield(optInfo,'offset'), optInfo.offset = []; end
    
    if ~isfield(optInfo,'score'), optInfo.score = []; end
    
    if isfield(fileInfo,'log') && ~isempty(fileInfo.log)
        
        logFile = fopen(fileInfo.log,'a');
        
        if logFile < 0, logFile = fopen('/dev/null'); end
        
    else
        
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
%{
    if ~isfield(optInfo,'initScore') || isempty(optInfo.initScore)
        
        disp('Calculating initial score...');
        jobPrefix = [fileInfo.prefix,'_init'];
        %simInfo.repeat = simInfo.predictRepeat;
        [jobList,dataList] = prepareJob(jobPrefix,[],[],optInfo.BS,...
            optInfo.CS,fileInfo,simInfo);
        optInfo.initScore = optInfo.funobj(jobList,dataList,paramFile);
        %cleanup(fileInfo,jobPrefix);
        disp('Initial score calculated.');
        
    end

    %if param_log == true

            optInfo.score = optInfo.initScore;
            optInfo.offset = zeros(1,nPara);

    %end
%}
    fprintf('Initial score: %e\n\n',optInfo.initScore);
    fprintf(logFile,'Initial score: %e\n\n',optInfo.initScore);
    
    %improvement = inf;
    s = optInfo.initGridSize;
    %nPara = size(optInfo.paraList,2);
    %quadFun = makeQuadFun(nPara);
    predictPara.BS = optInfo.BS;
    predictPara.CS = optInfo.CS;
    [x,xSize] = prepareGrid(s,nPara);
    y = zeros(xSize,2);
    y(:,2) = 3*optInfo.noise;
    minScore = optInfo.initScore;
    count = 0;
    noImprove = 0;
    
    while noImprove <= optInfo.maxIter %s >= optInfo.minGridSize
    %while r <= optInfo.maxBias % && improvement >= optInfo.eps
            
        if isfield(fileInfo,'debugFile1') && ...
                exist(fileInfo.debugFile1,'file')
            
            keyboard;
            
        end
        
        count = count + 1;
        %s = optInfo.getS(r);
        jobList = [];
        dataList = [];
        %simInfo.repeat = simInfo.sampleRepeat;
        iterPrefix = [fileInfo.prefix,'_iter',num2str(count)];
        fprintf('Search iteration: %d\nGrid size: %e\n',count,s);
        fprintf(logFile,'Search iteration: %d\nGrid size: %e\n',count,s);
        
        for i = 1 : xSize

            %if sum(abs(x(i,:))) == 0, continue; end

            fprintf('\rPreparing sample points: %d out of %d',i,xSize);
            [bs,cs] = setParameter(optInfo.BS,optInfo.CS,...
                optInfo.paraList,x(i,:));
            jobPrefix = [iterPrefix,'/grid',num2str(i)];
            [jobList,dataList] = prepareJob(jobPrefix,jobList,...
                dataList,bs,cs,fileInfo,simInfo);

        end
        
        fprintf('\n');
        %sendJob(jobList,qInfo);
        y(:,1) = optInfo.funobj(jobList,dataList,paramFile);
        
        %if param_log == true

            optInfo.score = [optInfo.score;y(:,1)];
            optInfo.offset = [optInfo.offset;x];
            save(fileInfo.paramFile,'fileInfo','qInfo','simInfo','optInfo');

        %end

        disp('Making prediction...');
        %offset = getOffset(x,[minScore;y],nPara,r,s,optInfo.lb,...
        %    optInfo.ub,quadFun);
        lb = ones(nPara,1) * optInfo.lb * s;
        ub = ones(nPara,1) * optInfo.ub * s;
        dx = (ub - lb) * optInfo.eps;
        params = struct('bounds',{lb,ub},'nreq',xSize,'p',0.5);
        
        if count == 1	

            [newPoints,xBest,yBest] = snobfit(tempFile,x,y,params,dx);
        
        else

            [newPoints,xBest,yBest] = snobfit(tempFile,x,y,params);

        end

        x = newPoints(:,1:nPara);
        cleanup(fileInfo,iterPrefix);
        %fprintf(logFile,'Sample scores:\n');
        %fprintf(logFile,'%e\t',[minScore;y]);
        fprintf(logFile,'Best offset of this iteration:\n');
        fprintf(logFile,'%e\t',xBest);
        fprintf('Best score of this iteration: %e\n',yBest);
        fprintf(logFile,'\nBest score of this iteration: %e\n',yBest);
%{        
        simInfo.repeat = simInfo.predictRepeat;
        [bs,cs] = setParameter(optInfo.BS,optInfo.CS,optInfo.paraList,...
            offset);
        jobPrefix = [fileInfo.prefix,'_predict',num2str(count)];
        [jobList,dataList] = prepareJob(jobPrefix,[],[],bs,cs,...
            fileInfo,simInfo);
        sendJob(jobList,qInfo);
        newScore = fitMultiCurve(dataList,simInfo);
        fprintf('Score of new prediction: %e\n',newScore);
        fprintf(logFile,'\nScore of new prediction: %e\n',newScore);
        %cleanup(fileInfo,jobPrefix);
%}        
        if yBest < minScore

            [bs,cs] = setParameter(optInfo.BS,optInfo.CS,...
                optInfo.paraList,xBest);
            predictPara.BS = bs;
            predictPara.CS = cs;
            %improvement = (minScore - newScore) / minScore;
            
            if minScore == inf 
                
                disp('Initialization finished.');
                fprintf(logFile,'Initial parameters:\n');
                
            else
                
                %s = s * 2;
                disp('Improved!');
                fprintf(logFile,'Improved! New parameters:\n');
                
            end
            
            for i = 1 : size(bs,1)
                
                fprintf(logFile,'%s %s %.15e %.15e\n',bs{i,1},bs{i,2},...
                    bs{i,3}(1),bs{i,3}(2));
                
            end
            
            for i = 1 : size(cs,1)
                
                fprintf(logFile,'%s %s %.15e\n',cs{i,1},cs{1,2},cs{i,3});
                
            end
            
            minScore = yBest;
            noImprove = 0;
                        
        else

            disp('No improvement. Decrease grid size...');
            %improvement = inf;
            %s = s / 2;
            noImprove = noImprove + 1;
            
        end
        
        fprintf(logFile,'\nNew requested offsets:\n');
        
        for i = 1 : xSize
            
            fprintf(logFile,'%.15e ',x(i,:));
            fprintf(logFile,'\n');
            
        end

        fprintf('Iteration %d finished.\n\n',count);
        fprintf(logFile,'Iteration %d finished.\n\n',count);
        
    end

    %if improvement < optInfo.eps
        
    %    disp('Minimal improvement reached.');
        
    %end
%{    
    if r > optInfo.maxBias
        
        disp('Minimal grid size reached.');
        
    end
%}
    fprintf('End search process for %s...\n\n',fileInfo.prefix);
    fprintf(logFile,'End search process for %s...\n\n',fileInfo.prefix);
    fclose(logFile);
    %predictPara.BS = optInfo.BS;
    %predictPara.CS = optInfo.CS;
    
end

function [ x, xSize ] = prepareGrid( s, nPara )

    if nPara == 1, x = [0;s;-s]; return, end
    
    xSize = 1+2*nPara^2;
    %x = zeros(1+2*nPara+4*nchoosek(nPara,2),nPara);
    %x = zeros(1+nPara+nPara^2,nPara);
    x = zeros(xSize,nPara);
    
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
%{
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
