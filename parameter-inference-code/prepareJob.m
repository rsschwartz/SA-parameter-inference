function [ jobList, dataList ] = prepareJob ( jobPrefix, jobList, dataList,...
    bindSite, confSwitch, fileInfo, simInfo)
%{
    if isfield(fileInfo,'preset'), preset = fileInfo.preset;
        
    else preset ='';
        
    end
%}    
    if isfield(simInfo,'maximer'), maximer = simInfo.maximer;
        
    else maximer = 200;
        
    end
    
    xmlTemplate = fileInfo.xmlRule;
    concNum = length(simInfo.conc);
    jobNum = simInfo.repeat*concNum;
    concFormat = ['%0',num2str(ceil(log(concNum+1)/log(10))),'d'];
    jobFormat = ['%0',num2str(ceil(log(simInfo.repeat+1)/log(10))),'d'];
    currJob = cell(1,jobNum);
    currData = cell(1,jobNum);
    xmlFile = cell(1,concNum);
    randList = unique(randi(100*jobNum,1,jobNum));
    
    while numel(randList) < jobNum
        
        randList = unique([randList,randi(100*jobNum,1,...
            jobNum-numel(randList))]);
        
    end
    
    bsConc = bindSite;
    bsNum = size(bindSite,1);
    
    for i = 1 : concNum
    
        for j = 1 : bsNum
            
            bsConc{j,3}(1) = bindSite{j,3}(1) + ...
                log(simInfo.conc(1)/simInfo.conc(i)) / log(10);
            
        end
        
        xmlFile{i} = [fileInfo.workFolder,jobPrefix,'/',num2str(i,...
            concFormat),'.xml'];
        modifyParameter(xmlTemplate,xmlFile{i},bsConc,confSwitch);
        
        for j = 1 : simInfo.repeat
        
            currPrefix = [jobPrefix,'/',num2str(j,jobFormat),'_conc',...
                num2str(i,concFormat)];
            currJob{(i-1)*simInfo.repeat+j} = [fileInfo.workFolder,...
                currPrefix,'.sh'];
            currData{(i-1)*simInfo.repeat+j} = [fileInfo.workFolder,...
                currPrefix,'.dat'];
            
            if isfield(fileInfo,'tempFolder')
                
                dataTemp = [fileInfo.tempFolder,...
                    strrep(jobPrefix,'/','_'),'_',num2str(j,jobFormat),...
                    '_conc',num2str(i,concFormat),'.dat'];
                fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');
                fprintf(fid,'%s %s %d %f %d %d > %s\n',fileInfo.command,...
                    xmlFile{i},simInfo.interval,simInfo.time,...
                    randList((i-1)*simInfo.repeat+j),maximer,dataTemp);
                fprintf(fid,'cp -bf %s %s\n',dataTemp,...
                    currData{(i-1)*simInfo.repeat+j});
                fprintf(fid,'cmp %s %s\nexport COPY_CMP=$?\n\n',...
                    dataTemp,currData{(i-1)*simInfo.repeat+j});
                fprintf(fid,'while [ $COPY_CMP -ne 0 ]\ndo\n');
                fprintf(fid,'cp -bf %s %s\n',dataTemp,...
                    currData{(i-1)*simInfo.repeat+j});
                fprintf(fid,'cmp %s %s\nexport COPY_CMP=$?\ndone\n\n',...
                    dataTemp,currData{(i-1)*simInfo.repeat+j});
                fprintf(fid,'rm -f %s\n',dataTemp);
                fclose(fid);
                
            else
                
                fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');
                fprintf(fid,'%s %s %d %f %d %d > %s\n',fileInfo.command,...
                    xmlFile{i},simInfo.interval,simInfo.time,...
                    randList((i-1)*simInfo.repeat+j),...
                    currData{(i-1)*simInfo.repeat+j});
                fclose(fid);
                
            end
            
        end
        
    end

    jobList = [jobList,currJob];
    dataList = [dataList,currData];
%{
    if isfield(fileInfo,'tempFolder')
        
        for i = 1 : concNum
            
            for j = 1 : simInfo.repeat
                
                dataTemp = [fileInfo.tempFolder,jobPrefix,'_',num2str(j,...
                    jobFormat),'_conc',num2str(i,concFormat),'.dat'];
                fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');
                fprintf(fid,'%s\n%s %s %d %f %d %d > %s\nmv %s %s\n',...
                    preset,fileInfo.command,xmlFile{i},simInfo.interval,...
                    simInfo.time,randList((i-1)*simInfo.repeat+j),...
                    maximer,dataTemp,dataTemp,currData{(i-1)*...
                    simInfo.repeat+j});
                fclose(fid);
                
            end
            
        end
    
    else
        
        for i = 1 : concNum
            
            for j = 1 : simInfo.repeat
                
                fid = fopen(currJob{(i-1)*simInfo.repeat+j},'w');
                fprintf(fid,'%s\n%s %s %d %f %d > %s\n',preset,...
                    fileInfo.command,xmlFile{i},simInfo.interval,...
                    simInfo.time,randList(i),currData{(i-1)*...
                    simInfo.repeat+j});
                fclose(fid);
                
            end
            
        end
        
    end
%}    
end

function [ xmlOut ] = modifyParameter ( xmlInFile, xmlOutFile, bindSiteList, confSwitchList )

    [xmlIn, xmlRoot] = xml_read(xmlInFile);
    pref.StructItem = false;
    
    if nargin < 3
        
        disp('Not enough parameters; nothing modified.');
        xml_write(xmlOutFile,xmlIn,xmlRoot,pref);
        
    else
    
        xmlOut = modifyBSTime(xmlIn,bindSiteList);
    
        if nargin >= 4
        
            xmlOut = modifyCSTime(xmlOut,confSwitchList);
        
        end

        xml_write(xmlOutFile,xmlOut,xmlRoot,pref);
        
    end
    
end

function [ xmlOut ] = modifyBSTime ( xmlIn, bindSiteList )

    xmlOut = xmlIn;
    bindingSiteType = xmlIn.BindingSiteTypes.BindingSiteType;
    numBindSite = size(bindingSiteType,1);
    numBSModify = size(bindSiteList,1);
    countModify = 0;
    
    for bsm = 1 : numBSModify
        
        bs1 = bindSiteList{bsm,1};
        bs2 = bindSiteList{bsm,2};
        timeModify = bindSiteList{bsm,3};
        numItemModify = length(timeModify);
        countModify = countModify + numItemModify;
        
        for bst = 1 : numBindSite
            
            currBS = bindingSiteType(bst).ATTRIBUTE.name;
            ptnNum = length(bindingSiteType(bst).Partner);
            
            for ptn = 1 : ptnNum
            
                partner = bindingSiteType(bst).Partner(ptn).ATTRIBUTE.name;
            
                if (strcmpi(currBS,bs1) && strcmpi(partner,bs2)) || ...
                        (strcmpi(currBS,bs2) && strcmpi(partner,bs1))
                
                    if numItemModify >= 1
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.bindTime = 10^timeModify(1);
                    
                    end
                
                    if numItemModify >= 2
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.breakTime = 10^timeModify(2);
                    
                    end
                
                    if numItemModify >= 3
                    
                        xmlOut.BindingSiteTypes.BindingSiteType(bst). ...
                            Partner(ptn).ATTRIBUTE.fastBindTime = 10^timeModify(3);
                    
                    end
                
                end
                
            end
            
        end
        
    end
%{
    if countModify == 0
            
        disp('No modification made on Bind/break/fastbind time.');
            
    else
            
        disp('Bind/break/fastbind time modified.');    
            
    end
%}    
end

function [ xmlOut ] = modifyCSTime ( xmlIn, confSwitchList )
    
    xmlOut = xmlIn;
    numConformation = size(xmlIn.ConformationalSwitch.ConformationTime,1);
    numCSModify = size(confSwitchList,1);
    
    for csm = 1 : numCSModify
        
        csFrom = confSwitchList{csm,1};
        csTo = confSwitchList{csm,2};
        timeModify = confSwitchList{csm,3};

        for c = 1 : numConformation
            
            currFrom = xmlIn.ConformationalSwitch.ConformationTime(c).ATTRIBUTE.name;
            
            if strcmpi(currFrom,csFrom) && length(timeModify) >= 1
                
                for cst = 1 : size(xmlIn.ConformationalSwitch.ConformationTime(c).List,1)
                    
                    currTo = xmlIn.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.name;
                    
                    if strcmpi(currTo,csTo)
                        
                        xmlOut.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.time ...
                            = 10^timeModify(1);
                        
                    end
                    
                end
                
            end
            
            if strcmpi(currFrom,csTo) && length(timeModify) == 2
                
                for cst = 1 : size(xmlIn.ConformationalSwitch.ConformationTime(c).List,1)
                    
                    currTo = xmlIn.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.name;
                    
                    if strcmpi(currTo,csFrom)
                        
                        xmlOut.ConformationalSwitch.ConformationTime(c).List(cst).ATTRIBUTE.time ...
                            = 10^timeModify(2);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
%{
    if numCSModify == 0
        
        disp('No modification made on Conformational Switch time.');
        
    else
        
        disp('Conformational Switch time modified.');
        
    end
%}    
end