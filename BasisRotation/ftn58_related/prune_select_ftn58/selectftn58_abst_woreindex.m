function reducedftn58=selectftn58_abst_woreindex(ftn58,orbitals,abst)
    % construct the ftn58 only with the orbitals 
    % the function will reindex the orbitals !
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    % the line will be dropped if hopping<abst
    
    reducedftn58=ftn58;
    
    % re-index
    [log2,ord2]=ismember(reducedftn58(:,2),orbitals);
    [log3,ord3]=ismember(reducedftn58(:,3),orbitals);
% I dont need reindex in kagome  
%     reducedftn58(:,2)=ord2;
%     reducedftn58(:,3)=ord3;
    
    % remove unwanted orbitals
    log = log2.*log3;
    oindices = find(~log);
    reducedftn58(oindices,:)=[];
    
    % remove hopping too small
    tindices=find(abs(reducedftn58(:,4))<abst);
    reducedftn58(tindices,:)=[];
    
    % compensate the first line in ftn58
    reducedftn58 = [[length(orbitals) size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return
