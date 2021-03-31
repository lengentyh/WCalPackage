function [reducedftn58,reindex_orbitals,ord2]=selectftn58(ftn58,orbitals)
    % construct the ftn58 only with the orbitals 
    % orbitals = [1:6,8,10]
    % Be attention to the order of the orbitals
    
    reducedftn58=ftn58;
    
    % re-index
    [log2,ord2]=ismember(reducedftn58(:,2),orbitals);
    [log3,ord3]=ismember(reducedftn58(:,3),orbitals);
    reducedftn58(:,2)=ord2;
    reducedftn58(:,3)=ord3;

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
    reducedftn58(indices,:)=[];
    
    % compensate the first line in ftn58
    reducedftn58 = [[size(orbitals,2) size(reducedftn58,1) 0 0 0 0 0]; reducedftn58]; 
    
return
