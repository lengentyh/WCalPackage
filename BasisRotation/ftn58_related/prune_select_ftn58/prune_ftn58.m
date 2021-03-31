function prunedftn58=prune_ftn58(ftn58,abst)
    % construct the ftn58 only with the orbitals 
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    % the line will be dropped if hopping<abst
    
    prunedftn58=ftn58;
    
    % remove hopping too small
    tindices=find(abs(prunedftn58(2:end,4))<abst);
    prunedftn58(tindices+1,:)=[];
    
    % update first line
    prunedftn58(1,2)=size(prunedftn58,1)-1;
    
return
