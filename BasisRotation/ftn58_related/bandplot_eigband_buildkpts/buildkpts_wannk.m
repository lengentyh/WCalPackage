function kpt=buildkpts(filename,subdirname)
%filename = 'YBCO123';
%subdirname = 'subdir10';

addpath(strcat('/home/lengentyh/WIEN2k/',filename,'/',subdirname));

b = strcat(subdirname,'_band.labelinfo.dat');
A = readwannk(b);
n=size(A.mat);
kpt = [];
i = 1;
for ord=1:A.mat(n(1,1),1)
    if ord > A.mat(i+1,1)
        i = i+1;
    end
    kpt(ord,1) = ord;
    kpt(ord,2) = ((A.mat(i+1,2) - A.mat(i,2))/(A.mat(i+1,1)-A.mat(i,1)))*(ord-A.mat(i,1)) + A.mat(i,2);
    kpt(ord,3) = ((A.mat(i+1,3) - A.mat(i,3))/(A.mat(i+1,1)-A.mat(i,1)))*(ord-A.mat(i,1)) + A.mat(i,3);
    kpt(ord,4) = ((A.mat(i+1,4) - A.mat(i,4))/(A.mat(i+1,1)-A.mat(i,1)))*(ord-A.mat(i,1)) + A.mat(i,4);
    kpt(ord,5) = ((A.mat(i+1,5) - A.mat(i,5))/(A.mat(i+1,1)-A.mat(i,1)))*(ord-A.mat(i,1)) + A.mat(i,5);
end

kpoints = kpt(:,3:5);
k4data.mat.kpoints = kpoints;
rmpath(strcat('/home/lengentyh/WIEN2k/',filename,'/',subdirname));
return