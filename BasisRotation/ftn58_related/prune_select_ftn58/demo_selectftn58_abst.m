load /home/lengentyh/Downloads/finerDOS/whole.mat ftn58
goaldir = strcat('/home/lengentyh/Downloads/finerDOS/','72x72_bulk_model');

planeorbitals = [4,5,6,7,8,9,10,11,12,13,19,20,21,22,23,24,25,26,27,28,29,30]; %[all-d all-p] 
planeftn58=selectftn58_abst(ftn58,planeorbitals,10^-5);
chainorbitals = [14,15,16,17,18,1,2,3,31,32,33,34,35,36]; % [all-d all-p vertical-p]
chainftn58=selectftn58_abst(ftn58,chainorbitals,10^-5);

cd(goaldir)
save 'ftn4prox_abst.mat' planeftn58 chainftn58;

allorbital = [planeorbitals, chainorbitals]; % make the orbitals ordered
realftn58 = selectftn58_abst(ftn58,allorbital,10^-5);
plnorbital = size(planeorbitals);
chnorbital = size(chainorbitals);

save 'ftn4prox_real_abst.mat' realftn58 plnorbital chnorbital;
