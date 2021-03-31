function reducedftn58=selectftn58_neighbor(ftn58,orbitals,n)
    % construct the ftn58 only with the neighboring orbitals 
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    % n = [2 2 1]
    % select only neighbors inside n cells (along both directions)
    
    reducedftn58=ftn58;
    
    % re-index
    [log2,ord2]=ismember(reducedftn58(:,2),orbitals);
    [log3,ord3]=ismember(reducedftn58(:,3),orbitals);
    reducedftn58(:,2)=ord2;
    reducedftn58(:,3)=ord3;
    
    % removes unwanted orbitals
    log = log2.*log3;
    oindices = find(~log);
    reducedftn58(oindices,:)=[];

    % removes neighbors too far
    neix = [-n(1):n(1)];  
    neiy = [-n(2):n(2)]; 
    neiz = [-n(3):n(3)]; 
    [nxlog,nxord] = ismember(reducedftn58(:,5),neix);
    [nylog,nyord] = ismember(reducedftn58(:,6),neiy);
    [nzlog,nzord] = ismember(reducedftn58(:,7),neiz);
    nlog = nxlog.*nylog.*nzlog;
    nindices = find(~nlog);
    reducedftn58(nindices,:)=[];

    % compensate the first line in ftn58
    reducedftn58 = [[size(orbitals,2) size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return
