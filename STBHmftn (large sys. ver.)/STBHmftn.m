%%%            Program for Constructing Slab ftn58sparse             %%%
%%% ---------------------------------------------------------------- %%%
%%% Necessary Input files:                                           %%%
%%% 1) ftn58sparse.mat                                               %%%
%%% 2) slab_info.mat (generated by SlabCell.m)                       %%%
%%% ---------------------------------------------------------------- %%%
%%% Input data description:                                          %%%
%%% nlayer    ==> number of layers in the slab supercell             %%%
%%% isCutSlab ==> if true, cutting a chosen surface                  %%%
%%%          ("uplimit" and "botlimit" are set for indicating your   %%%
%%%           up and bottom coordinates along hkl direction;         %%%
%%%           "xyz_dir" indicates the xyz cartesian axis)            %%%
%%% ---------------------------------------------------------------- %%%
% clear all
function STBHmftn()
wcal = ReadInput('input.txt');

%%% Input this manually %%%
isRT = 0;
RT   = [1 -1 0;1 1 0;0 0 sqrt(2)]/sqrt(2);
load(wcal.ref);
load('slab_info.mat');
%%% ------------------- %%%

tic
%% --- Read Input information --- %%%
dd  = ftn58sparse.dd*bulkbasis*transpose(inv(transform_matrix));
ib1 = ftn58sparse.ij(:,1);
ib2 = ftn58sparse.ij(:,2);
tt  = ftn58sparse.tt;

ijs =zeros(size(ib1,1)*500,2);
tts =zeros(size(ib1,1)*500,1);
dds =zeros(size(ib1,1)*500,3);

%% --- Actual Procedure --- %%%
fid = fopen('slab_ftn58.dat','w');
kk=1;
for ii=1:size(orbitps, 1) % total number of needed orbitals (in a unit cell) 
    if mod(ii,50)==0
        fprintf('# of Orbital = %4i/%i\n',ii,size(orbitps,1));
    end
    temp1 = find(ib1==orbitps(ii,2)); %look for all possible translations of obrital(i)
    for jj=1:size(temp1, 1) % run over all B
        temp2 = find(orbitps(:,2)==ib2(temp1(jj))); %look for hopping parameters of orbital(i) and all required orbitals ###
        if(temp2>0)
            component = orbitps(ii, 4:6) + orbitps(temp2(1),7:9) - orbitps(ii,7:9) + dd(temp1(jj),:);
                   % hopping_AB = orb(new_bulk frac.)_A + atom(new_bulk frac.)_B -
                   % atom(new_bulk frac.)_A + ori_translation(new_bulk frac.)
            index = find((abs(orbitps(1:end,6)-component(3))<0.001)&(orbitps(:,2)==ib2(temp1(jj))));
                % search the orb that corresponding to B   
            if(index>0)
                %ijs(kk,1:2) = [ii, orbitps(index,1)];
                %tts(kk)     = tt(temp1(jj));
                %dds(kk,:)   = round(component-orbitps(index, 4:6));
                % hopping_frac_AB = hopping_AB - Dis_B
				
				% using is slow for large array
				fprintf(fid,"%5d\t%5d\t%10.8f\t%10.8f\t%5d\t%5d\t%5d\n",ii,orbitps(index,1),real(tt(temp1(jj))),imag(tt(temp1(jj))),round(component-orbitps(index,4:6)));
                kk = kk+1;
            end
        end
    end
end

disp(kk)
fclose(fid);
toc

%% --- reread the ij, tt, dd --- %%%
ftn58_temp = [];

%j=0;
fid        		= fopen('slab_ftn58.dat','r');
format     		= '%d %d %f %f %d %d %d';
sizeftn58  		= [7 Inf];
ftn58      		= fscanf(fid,format,sizeftn58);
ftn58_temp      = ftn58';
ftn58_temp(:,3) = ftn58_temp(:,3) + ftn58_temp(:,4);
ftn58_temp(:,4) = [];
%while ~feof(fid)
%    j=j+1;
%    tline=fgetl(fid);
%
%        i=1;
%        while tline(i)==' ' %skip blank
%            i=i+1;
%        end
%
%    temp       = split(tline);
%    ij_temp    = str2double(temp(1:2,1))';
%    tt_temp    = str2double(temp(3,1)) + 1i*str2double(temp(4,1));
%    dd_temp    = str2double(temp(5:7,1))';
%    ftn58_temp(j,:) = [ij_temp tt_temp dd_temp];
%end
fclose(fid);

ijs = ftn58_temp(:,1:2);
tts = ftn58_temp(:,3);
dds = ftn58_temp(:,4:6);

%% --- Rotating the new unit vectors --- %%
if (isRT)
    BR            = BR*RT'; 
end

%% --- Save as Sftn58sparse format --- %%%
for i=1:ftn58sparse.Nat
    aname = getfield(ftn58sparse.Ainfo, {i}, 'Atom');
end

for i=1:size(superatomps,1)
    atominfo(i).Atom(1,:) = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'Atom');
    atominfo(i).Position(1,1:3) = superatomps(i,9:11);
    atominfo(i).Norb = superatomps(i,4);
    atominfo(i).OrbitIndex = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'OrbitIndex');
    atominfo(i).Orbit(1,:) = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'Orbit');
    atominfo(i).OrbitID    = getfield(ftn58sparse.Ainfo, {superatomps(i,5)}, 'OrbitID');
end

% ob1 = [1:6 7:12 19:24 25:30];
% ob2 = [7:12 13:18 25:30 31:36];
% ttR = 2.5;
% d00 = find(dds(:,1)==0&dds(:,2)==0);
% ii  = ijs(:,1);
% jj  = ijs(:,2);
% for i=1:length(ob1)/3
%     id_ob1{i}  = ob1((i-1)*3+1:i*3);
%     id_ob2{i}  = ob2((i-1)*3+1:i*3);
% end
% for i=1:length(id_ob1)
%     for ib1=1:length(id_ob1{1,1})
%         for ib2=1:length(id_ob2{1,1})
%             hop1 = find(ii(d00)==id_ob1{1,i}(ib1)&jj(d00)==id_ob2{1,i}(ib2));
%             hop2 = find(ii(d00)==id_ob1{1,i}(ib2)&jj(d00)==id_ob2{1,i}(ib1)); 
%             tts(d00(hop1)) = tts(d00(hop1))*ttR;
%             tts(d00(hop2)) = tts(d00(hop2))*ttR;
%         end
%     end
% end

Sftn58sparse.System  = 'Slab_ftn58';
Sftn58sparse.Ainfo   = atominfo;
Sftn58sparse.hkl     = hkl;
Sftn58sparse.abc     = ftn58sparse.abc;
Sftn58sparse.BR3D    = ftn58sparse.BR;
Sftn58sparse.BR2D    = BR;
Sftn58sparse.Nat     = length(atominfo); 
Sftn58sparse.ver     = 'type1';
Sftn58sparse.isSO    = ftn58sparse.isSO;
Sftn58sparse.norb    = size(orbitps,1);
Sftn58sparse.Orbitps = [orbitps(1:end,1:3) orbitps(1:end,4:6) orbitps(1:end,10)];
Sftn58sparse.ij      = ijs;
Sftn58sparse.tt      = tts;
Sftn58sparse.dd      = dds(:,1:2);

save Sftn58sparse.mat Sftn58sparse -v7.3
