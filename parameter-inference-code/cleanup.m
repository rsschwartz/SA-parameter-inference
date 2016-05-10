function cleanup ( fileInfo, prefix )

    errFile = [];
    
    if prefix(end) == '/', prefix = prefix(1:end-1); end
    
    if isfield(fileInfo,'errFolder') && ...
            exist(fileInfo.errFolder,'dir') == 7
    
        %[status,~] = system(['cd ',fileInfo.errFolder,';ls ',prefix,...
        %    '_*']);
        
        if fileInfo.errFolder(end) ~= '/'
            
            fileInfo.errFolder = [fileInfo.errFolder,'/'];
            
        end
        
        %[status,~] = system(['cd ',fileInfo.errFolder,';ls *']);
        
        [status,~] = system(['ls ',fileInfo.errFolder,'*']);
                
        if status == 0
        
            errFile = ['err_',prefix,'.tgz'];
            %[status,msg] = system(['cd ',fileInfo.errFolder,...
            %    ';tar -czf ',fileInfo.workFolder,errFile,' ',prefix,'_*']);
            
            [status,msg] = system(['cd ',fileInfo.errFolder,...
                ';tar -czf ',fileInfo.workFolder,errFile,' *']);
    
            if status == 0
        
                system(['rm -rf ',fileInfo.errFolder,'*']);
        
            else
        
                disp('Error creating tar ball!')
                disp(errFile);
                disp(msg);
                
            end
        
        end
        
    end

    tgzFile = [fileInfo.tempFolder,prefix,'.tgz'];
    [status,msg] = system(['cd ',fileInfo.workFolder,';tar -czf ',...
        tgzFile,' ',prefix,'* ',errFile]);
    
    if status == 0
        
        system(['rm -rf ',fileInfo.workFolder,prefix,'*']);
        system(['cp -bf ',tgzFile,' ',fileInfo.workFolder]);
        status = system(['cmp ',tgzFile,' ',fileInfo.workFolder,prefix,...
            '.tgz']);
        
        while status ~= 0
            
            system(['cp -bf ',tgzFile,' ',fileInfo.workFolder]);
            status = system(['cmp ',tgzFile,' ',fileInfo.workFolder,...
                prefix,'.tgz']);
            
        end
        
        system(['rm -f ',tgzFile]);
        
        if ~isempty(errFile)
            
            system(['rm -f ',fileInfo.workFolder,errFile]);
            
        end
        
    else
        
        disp('Error creating tar ball!')
        disp(tgzFile);
        disp(msg);
        
    end
    
end