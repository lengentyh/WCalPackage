function reducedftn58=selectftn58_left_right(ftn58,lorbitals,rorbitals)
    % construct the ftn58 only with the orbitals 
    % the orbitals are reindex
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    
    reducedftn58=ftn58;
    
    % re-index
    [log2,ord2]=ismember(reducedftn58(:,2),lorbitals);
    [log3,ord3]=ismember(reducedftn58(:,3),rorbitals);
    reducedftn58(:,2)=ord2;
    reducedftn58(:,3)=ord3;
    
    % remove unwanted orbitals
    log = log2.*log3;
    indices = find(~log);
    reducedftn58(indices,:)=[];
    
    % compensate the first line in ftn58
%     reducedftn58 = [[2 size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return