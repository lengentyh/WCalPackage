function [Ek,p,PCD] = bandplot_concise_ftn58sparse_ex(ftn58sparse, Ef, isSP, kpt_filename, PCD_orb)
% isSP = 0 | 1;
% kpt_filename : the name of the file '.labelinfo.dat'
% more detail information please see the help of 'readkpt.m'
% PCD_orb : orbitals that you want to calculate their partial charge distributions
% e.g. PCD_orb = [1 2 3];
% PCD_orb = [] to NOT have the calculation about partial charge distributions
% add PCD_orb_2 to plot another partial charge distribution in case ... etc.

if isSP ~= 0 & isSP ~= 1
	fprintf('Please gives correct information about spin polarization ! \n');
	return
end 

%% Open a pool for parallel computation %%%
% --------------------------------------- %
delete(gcp);
pool=parpool('local');
pool.IdleTimeout=1800; %minutes
%=======================================================================================================
%% ----------------- Initialization --------------------- %%%
% --------------------------------------------------------- %
norb = ftn58sparse.norb;
ii   = ftn58sparse.ij(:,1);
jj   = ftn58sparse.ij(:,2);
dd   = ftn58sparse.dd;
tt   = ftn58sparse.tt;
BR   = ftn58sparse.BR;
abc  = ftn58sparse.abc;
Sz   = [1 0;0 -1];
Sy   = [0 -1i;1i 0];
Sx   = [0 1;1 0];
%=======================================================================================================
%% --- Read The K-path & Generate the Plotting Data ----- %%%
% --------------------------------------------------------- %
b = strcat(kpt_filename,'.labelinfo.dat');
A = readkpt(b);
%------------------------------------------
% A.label(j)= label;
% A.mat(j,:)= [partition;mesh;position];
%------------------------------------------
n = size(A.mat,1);
symlb = A.label;
for j = 1:n
	if j<= n-1
		p_temp   = linspace(A.mat(j,5),A.mat(j+1,5),A.mat(j+1,1));
		kpt_temp = [linspace(A.mat(j,2),A.mat(j+1,2),A.mat(j+1,1))' ...
					linspace(A.mat(j,3),A.mat(j+1,3),A.mat(j+1,1))' ...
					linspace(A.mat(j,4),A.mat(j+1,4),A.mat(j+1,1))' ];
	end				
	if j == 1
		p     = p_temp;
		kpt   = kpt_temp;
		sympt = A.mat(j,5);
	else
		p     = [p(1:end-1),p_temp];
		kpt   = [kpt(1:end-1,:);kpt_temp];
		sympt = [sympt, A.mat(j,5)];
	end
end
% --------------------------------------------------------- %
%% --- Paragraph Output: 'p', 'kpt', 'sympt', 'symlb' --- %%%
% --------------------------------------------------------- %
% p     : the list records the position of each k-points on the axis of the dispersion plot
% kpt   : the list records the k-points along the choosen k-path
% sympt : the list records the position of the high symmetry points on the axis of the dispersion plot
% symlb : the list records the corresponding characters of the high symmetry points
% ------------------------------------------------------------------------------------------------------
% p     = [ p1, p2, ... ];
% kpt   = [ [kpt1.x kpt1.y kpt1.z]; [kpt2.x kpt2.y kpt2.z]; ... ];
% sympt = [ sympt1, sympt2, ... ];
% symlb = [ 'Gamma', 'X', ... ];
%=======================================================================================================
%% --- Calculate the Eigenvalues, Eigenvectors and the Partial Charge Distriution ---- %%%
% -------------------------------------------------------------------------------------- %
% ---------------------- Preallocation ---------------------- % 
% ----------------------------------------------------------- % 
nks     = size(kpt,1);
kpoints = 2*pi*kpt;
Ek      = zeros(nks,norb);
PCD 	= zeros(nks,norb);
PCD_exp = zeros(norb,norb);
% concepts to plot another partial charge distribution
% PCD_2     = zeros(nks,norb);
% PCD_exp_2 = zeros(nks,norb);
% ----------------------------------------------------------- %
% --- Expectation Value for the partial charge distribution - % 
% --- Using block diagonalized Identity matrix -------------- %
% ----------------------------------------------------------- %
if length(PCD_orb) ~= 0
	v              = ones(length(PCD_orb),1);
	PCD_exp_sparse = sparse(PCD_orb,PCD_orb,v,norb,norb);
	PCD_exp        = full(PCD_exp_sparse); 
end
% ----------------------------------------------------------- %
% ---------------------- Calculation ------------------------ % 
% ----------------------------------------------------------- %
tic
fprintf('Start calculating ... \n')
parfor ik=1:nks
    Hsparse   = sparse(ii,jj,exp(1i*dd*kpoints(ik,:)').*tt,norb,norb);
    HH        = full(Hsparse);
    HH        = (HH+HH')/2;
    [vec, Etemp] = eig(HH);
     
	% Record the calculation 
    % Ham{ik,1}    = HH;
    % eigvec{ik,1} = vec;
    
    Ek(ik,:)  = diag(Etemp);
	% Calculate the partial charge distributions
    if length(PCD_orb) ~= 0
		PCD(ik,:) = diag(vec' * PCD_exp * vec);
    end 

end
Ek = Ek - Ef;
toc
%=======================================================================================================
%close the pool
delete(gcp);

return

