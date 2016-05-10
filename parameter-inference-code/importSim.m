function [ simData ] = importSim ( simFile )

    simData = [];
    
    try
        
        fid = fopen(simFile,'r');
        
    catch fopenErr
        
        fprintf('Failed to open file: %s\n',simFile);
        throw(fopenErr);
        
    end
    
    if fid < 0, fprintf('Failed to open file: %s\n',simFile); return, end
        
    count = 0;
    str = fgetl(fid);
    
    if isempty(str) || str(1) == -1
        
        fprintf('Error opening file: %s\n',simFile);
        
        return
    
    end
    
    while ~isempty(str) && ~isNum(str(1)) && ~feof(fid)
        
        str = fgetl(fid);
        count = count + 1;
        
    end
    
    if feof(fid), return, end
    
    fclose(fid);
    simData = importdata(simFile,' ',count);
    
    if isfield(simData,'data'), simData = simData.data; end

end

function [ bool ] = isNum ( aChar )

    if aChar >= 48 && aChar <=57, bool = true;
    
    else bool = false;
        
    end

end