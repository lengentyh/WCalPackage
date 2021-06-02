function new_ftn58sparse = basis_transfromation_ftn58sparse(ftn58sparse,V)
% Transform the basis of ftn58saprse according to V
% The ftn58spare must in Bloch representation


% initialization 
new_ftn58sparse = ftn58sparse;
ij		= ftn58sparse.ij;
ii  		= ftn58sparse.ij(:,1);
jj		= ftn58sparse.ij(:,2);
tt   		= ftn58sparse.tt;
dd   		= ftn58sparse.dd;
norb            = ftn58sparse.norb;
nbond           = size(ij,1);

% dummy evasion
if length(find(mod(ftn58sparse.dd,1))) ~= 0
	fprintf('Please transform the ftn58sparse into Bloch representation ! \n')
end
if size(V,1) ~= norb | size(V,2) ~= norb
	fprintf('Please check the size of the transformation matrix ! \n')
	return
end
if abs(V'*V - eye(norb)) > 1.0e-10
	fprintf('Please check the unitarity of the transformation matrix ! \n')
	return
end

% add column in ftn58 to classify translation
trsl		= unique(dd,'rows');
[logn,ordn] 	= ismember(dd,trsl,'rows');
    
% pull out and construct hopping matrix
nlen 		= length(trsl);
t 		= zeros(norb,norb,nlen);
for j=1:nbond  
    	m=ij(j,1);
    	k=ij(j,2);
	n=ordn(j);
	t(m,k,n)=t(m,k,n)+tt(j);
end

% transform (t) for each cell
trot = zeros(norb,norb,nlen);
ROT  = V;

ij_temp = zeros(1,2);
tt_temp = zeros(1);
dd_temp = zeros(1,3);
for j=1:nlen
	trot(:,:,j) = ROT'*t(:,:,j)*ROT;
	
	% put back trot into ftn58sparse
	trot_sparse = sparse(trot(:,:,j));
	[rol,col]   = find(trot_sparse);
	
	ij_temp = [ij_temp;[rol,col]];
	tt_temp = [tt_temp;nonzeros(trot_sparse)];
	dd_temp = [dd_temp;repmat(trsl(j,:),length(rol),1)];
end
new_ftn58sparse.ij = ij_temp(2:end,:);
new_ftn58sparse.tt = tt_temp(2:end);
new_ftn58sparse.dd = dd_temp(2:end,:);

