function [ bs, cs ] = setParameter ( bsRef, csRef, paraList, x )

    if exist('bsRef','var') && ~isempty(bsRef)
    
        bsTime = cell2mat(bsRef(:,3));
        bs = bsRef;
        
    else bs = [];
        
    end
    
    if exist('csRef','var') && ~isempty(csRef)
        
        csTime = cell2mat(csRef(:,3));
        cs = csRef;
        
    else cs = [];
        
    end
    
    for i = 1 : size(paraList,2);
        
        if paraList{i}{1} == 'b' && ~isempty(bs)
            
            bsPara1 = paraList{i}{2};
            bsPara2 = paraList{i}{3};
            bsTime(bsPara1,bsPara2) = bsTime(bsPara1,bsPara2) + x(i);
            
        end

        if paraList{i}{1} == 'c' && ~isempty(cs)
            
            csPara = paraList{i}{2};
            %csPara2 = paraList{i}{3};
            csTime(csPara) = csTime(csPara) + x(i);
            
        end

    end

    if ~isempty(bs)

        for i = 1 : size(bsTime,1), bs{i,3} = bsTime(i,:); end
        
    end
    
    if ~isempty(cs)

        for i = 1 : size(csTime,1), cs{i,3} = csTime(i,:); end
        
    end
    
end