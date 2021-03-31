load /home/lengentyh/Downloads/kagome/ftn58sparse.mat ftn58sparse
ftn58 = ftn58sparse2ftn58(ftn58sparse);
abst=1e-3;
ftn58 = prune_ftn58(ftn58,abst);

upftn58 = ftn58;
% use only spin up bands
log1 = find(upftn58(:,2)>(upftn58(1,1)/2));
upftn58(log1,:)=[];
log2 = find(upftn58(:,3)>(upftn58(1,1)/2));
upftn58(log2,:)=[];

dnftn58 = ftn58;
% use only spin up bands
log3 = find(dnftn58(:,2)<=(dnftn58(1,1)/2));
dnftn58(log3,:)=[];
log4 = find(dnftn58(:,3)<=(dnftn58(1,1)/2));
dnftn58(log4,:)=[];

save('spinftn58.mat','upftn58','dnftn58')