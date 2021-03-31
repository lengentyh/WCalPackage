function reducedftn58=selectftn58_nondiag_abst(ftn58,orbitals1,orbitals2,abst)
    % construct the non-diagonal part (orbitals1, orbitals2) of ftn58 
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    % the line will be dropped if hopping<abst
    
    reducedftn58=ftn58;
    
    % compensate the lower triangle
    reducedftn58 = ftn58;
    for j=2:size(ftn58,1)
        if ftn58(j,2)<ftn58(j,3)
            reducedftn58 = [reducedftn58; [1+size(reducedftn58,1) ftn58(j,3) ftn58(j,2) conj(ftn58(j,4)) ftn58(j,5:7)] ];
        end
    end
    
    % re-index
    [log2,ord2]=ismember(reducedftn58(:,2),orbitals1);
    [log3,ord3]=ismember(reducedftn58(:,3),orbitals2);
    reducedftn58(:,2)=ord2;
    reducedftn58(:,3)=ord3;
    
    % remove unwanted orbitals
    log = log2.*log3;
    oindices = find(~log);
    reducedftn58(oindices,:)=[];
    
    % remove hopping too small
    tindices=find(abs(reducedftn58(:,4))<abst);
    reducedftn58(tindices,:)=[];
    
    % compensate the first line in ftn58
    reducedftn58 = [[size(orbitals1)+size(orbitals2) size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return