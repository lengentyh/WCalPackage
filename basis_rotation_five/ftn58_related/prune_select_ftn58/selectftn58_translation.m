function reducedftn58=selectftn58_translation(ftn58,n)
    % select only nth cell in ftn58
    % n = [1x3]
    
    reducedftn58=ftn58;
    
    % re-index
    [log,ord]=ismember([reducedftn58(:,5:7)],n,'row');
    
    % removes unwanted orbitals
    oindices = find(~log);
    reducedftn58(oindices,:)=[];

    % compensate the first line in ftn58
    reducedftn58 = [[ftn58(1,1) size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return
