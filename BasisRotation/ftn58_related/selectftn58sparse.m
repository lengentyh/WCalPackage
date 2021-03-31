function [reducedftn58sparse,reindex_orbitals,ord2]=selectftn58sparse(ftn58sparse,orbitals)
    % construct the ftn58 only with the orbitals 
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    
    % initialization
    ij = ftn58sparse.ij;
    tt = ftn58sparse.tt;
    dd = ftn58sparse.dd;
    reducedftn58sparse = ftn58sparse;
    
    % re-index
    [log2,ord2]=ismember(ij(:,1),orbitals);
    [log3,ord3]=ismember(ij(:,2),orbitals);
    ij(:,1)=ord2;
    ij(:,2)=ord3;

%    % debug
%    c = 1;
%    for ord = 1:length(ord3)
%        if ord3(ord) ~= 0
%    		reindex_orbitals(c)=orbitals(ord3(ord));
%		c = c + 1;
%	end
%    end
    
    % remove unwanted orbitals
    log = log2.*log3;
%     for j=1:size(reducedftn58,1)
%         if log(j)==0
%             reducedftn58(j,:)=[];
%         end
%     end
    indices = find(~log);
    ij(indices,:) = [];
    dd(indices,:) = [];
    tt(indices)   = [];

    reducedftn58sparse.ij = ij;
    reducedftn58sparse.dd = dd;
    reducedftn58sparse.tt = tt;
    reducedftn58sparse.norb = length(orbitals);
    
return
