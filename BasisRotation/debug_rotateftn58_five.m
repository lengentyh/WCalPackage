clear all
addpath ./ftn58_related
addpath ./bandplot_concise

% use the TB model of Tb(MnSn)6 
load ftn58sparse.mat
norb = ftn58sparse_temp.norb;

% ====================================================================================================================
% INPUTS
orbital 	= 42:44;  % the orbitals to be rotated
diff_tol_Ek     = 10^-14; % the tolerance for difference in band 			( before and after rotation )
diff_tol_PCD    = 10^-12; % the tolerance for difference in Partial Charge Distribution ( before and after rotation )
PCD_orb_1 	= [43];   % e.g. Sn, px
PCD_orb_2 	= [44];	  % e.g. Sn, py
Ef      	= 6.1679; % the Fermi level
% =====================================================================================================================

% selsect certain translations
[Lia,Lob] = ismember(ftn58sparse.dd,[-5 -3 -2],'rows');
imag_Lia = find(~Lia);
ftn58sparse_temp = ftn58sparse;
ftn58sparse_temp.dd(imag_Lia,:) = [];
ftn58sparse_temp.tt(imag_Lia,:) = [];
ftn58sparse_temp.ij(imag_Lia,:) = [];

% select the spin up part
select_up_orbital = 1:59;
ftn58sparse_temp = selectftn58sparse(ftn58sparse,select_up_orbital);
fprintf('selection is over \n')

% ================================================================================
% ROTATE the ftn58sparse
[new_ftn58sparse,trot] = rotateftn58_five(ftn58sparse_temp,pi/2,0,0,'p',orbital);
fprintf('rotation is over \n')
% ================================================================================

% calculate the band structure to check
% that the operation is indeed a rotation of bais
tic
[Ek_origin,p,PCD_origin_1] = bandplot_concise_ftn58sparse_ex(ftn58sparse_temp, Ef, 0, 'kpt', PCD_orb_1);
[Ek_rot,~,PCD_rot_1] = bandplot_concise_ftn58sparse_ex(new_ftn58sparse, Ef, 0, 'kpt', PCD_orb_1);
[~,~,PCD_origin_2] = bandplot_concise_ftn58sparse_ex(ftn58sparse_temp, Ef, 0, 'kpt', PCD_orb_2);
[~,~,PCD_rot_2] = bandplot_concise_ftn58sparse_ex(new_ftn58sparse, Ef, 0, 'kpt', PCD_orb_2);
toc
fprintf('Ek calculations are over \n')

% visualize the results
figure('Name','partial charge distribution org')
hold on
tic
plot(p,Ek_origin,'k')
for j = 1:norb
	if length(PCD_orb_1) ~= 0 && length(PCD_orb_2) ~= 0
		scatter(p,Ek_origin(:,j),100*PCD_origin_1(:,j)+1.0e-10,'r','filled');
		scatter(p,Ek_origin(:,j),100*PCD_origin_2(:,j)+1.0e-10,'b','filled');
	end
	drawnow
end
hold off
title('original')
set(gcf, 'Position',  [150, 150, 2000, 1600])
% ------------------------------------------------------------------------------------
figure('Name','partial charge distribution rot')
hold on
plot(p,Ek_rot,'k')
for j = 1:norb
	if length(PCD_orb_1) ~= 0 && length(PCD_orb_2) ~= 0
		scatter(p,Ek_rot(:,j),100*PCD_rot_1(:,j)+1.0e-10,'r','filled');
		scatter(p,Ek_rot(:,j),100*PCD_rot_2(:,j)+1.0e-10,'b','filled');
	end
	drawnow
end
hold off
title('rotated')
set(gcf, 'Position',  [150, 150, 2000, 1600])
drawnow
toc
fprintf('plots are over \n')

% ----------------------------------------------------------------------
% PRINT the QUANTITIVE results for checking band structure
Lia_Ek        = ismembertol(Ek_rot,Ek_origin,diff_tol_Ek);
bad_pts_Ek    = find(~Lia_Ek);
bad_pts_PCD_1  = find(abs(PCD_origin_1 - PCD_rot_2)>diff_tol_PCD);
bad_pts_PCD_2  = find(abs(PCD_origin_2 - PCD_rot_1)>diff_tol_PCD);
fprintf('Amount of bad points of Ek %f \n',length(bad_pts_Ek))
fprintf('Amount of bad points of PCD_orb_1 %f \n',length(bad_pts_PCD_1))
fprintf('Amount of bad points of PCD_orb_2 %f \n',length(bad_pts_PCD_2))

% ======================================================================
% THE FATAL TEST of the rotation
% check the hermiticity of the ftn58sparse
new_ftn58 = ftn58sparse2ftn58(new_ftn58sparse);
check_new = check_hopping_diagonal(new_ftn58,diff_tol_Ek);
fprintf('the check is over \n')
% =======================================================================

rmpath ./ftn58_related
rmpath ./bandplot_concise
