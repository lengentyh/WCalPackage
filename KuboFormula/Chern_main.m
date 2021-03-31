%%% Program for calculating Spin Hall conductivity via Kubo formula (3D) %%%
%%% -------------------------------------------------------------------- %%%
%%% Define spin Berry curvature Omega^k_(i,j), where i and j stands for  %%%
%%% the direction of spin current and electric field. And k is for spin  %%%
%%% direction (k=x,y,z).                                                 %%%
%%% Example: Sig^{Sk}_{Ji,Ej}                                            %%%
%%% -------------------------------------------------------------------- %%%
clear all

%% Initial setup for i, j, and k  %%%
Nk    = [50,50];    %--- k-mesh
isRT  = 0;          %--- is the system beening Rotation-Transformed or not
Ktype = 1;          %--- is the mesh along cover part of BZ or along a path

Temp  = 1e-3;
Ef    = -3.7016;
Emu   = [-0.1,0.1];
nE    = 20; 

load ftn58sparse_soc.mat
norb = ftn58sparse.norb;
%% K-mesh setup (at kz == 0 as the enclosed surface)%%%
if Ktype==1
    N1 = Nk(1); N2 = Nk(2);
    p1 = linspace(0,1,N1+1); 
    p2 = linspace(0,1,N2+1); 
    p1 = p1(1:N1); 
    p2 = p2(1:N2); 
    [K1,K2] = meshgrid(p1,p2);
    kpts = zeros(N1*N2,3);

%     BR  = ftn58sparse.BR;
    ik = 1;
    for ii=1:N1
        for jj=1:N2
                kpts(ik,:) = [K1(ii,jj) K2(ii,jj) 0]*2*pi ;
                ik = ik + 1;
        end
    end
    nks = length(kpts);
end

%% K-path along high symmetry path %%%
if Ktype==0
    nk = Nk(1);
    p1 = [linspace(2/3,1/3,nk+1)' linspace(-1/3,1/3,nk+1)' linspace(0.0,0.0,nk+1)'];
    %p1 = [linspace(0.0,0.0,nk+1)' linspace(0.0,0.5,nk+1)' linspace(0.0,0.0,nk+1)'];
    %p2 = [linspace(0.0,1/3,nk)' linspace(0.5,1/3,nk)' linspace(0.0,0.0,nk)'];
    %p3 = [linspace(1/3,2/3,nk)' linspace(1/3,-1/3,nk)' linspace(0.0,0.0,nk)'];
    %p4 = [linspace(2/3,0.0,nk)' linspace(-1/3,0.0,nk)' linspace(0.0,0.0,nk)'];
    %kpts = [p1(1:nk,:);p2(1:nk,:);p3(1:nk,:);p4(1:nk,:);]*2*pi;
    kpts = [p1(1:nk,:)]*2*pi;
    nks  = length(kpts);
end

%% --- Berry phase calculation --- %%%
tic
if isRT ==0
    if Ktype==1
        Berry = zeros(nks,norb,9);
        Ek    = zeros(nks,norb);
        tic
        for nk=1:nks
            fprintf('total k=%d, now=%d\n',nks,nk)
%             time_init=tic;
%             for n2=1:N2
%                 kpoints(1:3)   = kpts(nk,:);
                [Btmp,Etmp]    = Chernxyz(kpts(nk,:),ftn58sparse);
                Berry(nk,:,:)  = [Btmp];
                Ek(nk,:)       = Etmp;               
%             end
%             if mod(nk,100)==0
%                 fprintf('%3i/%i: %.3fs [%.2fm]\n',n1,N1,toc,toc(time_init)/60);
%             end
        end    
        save Berry.mat Berry Ek
    else
        Berry = zeros(nks,norb,9);
        Ek    = zeros(nks,norb);
        tic
       for ik=1:nks
            disp(ik);
%             kpoints       = kpts(ik,:);
            [Btmp, Etmp]  = Chernxyz(kpts(ik,:),ftn58sparse);
            Berry(ik,:,:) = [Btmp];
            Ek(ik,:)      = Etmp;
        end
        save Berry_kpath.mat Berry Ek
    end
else
    load Berry_mx.mat
end

toc
return
%% --- QSHE at various energies --- %%%
load Berry.mat
cquantum = 1; % the unit becomes (e^2/h)  
kB       = 8.6173324E-5;
BR       = ftn58sparse.BR;
abc      = ftn58sparse.abc;

if Ktype==1
    T       = [BR(:,1)*abc(1) BR(:,2)*abc(2) BR(:,3)*abc(3)];
    %Volume  = abs(T(3,:)*cross(T(1,:),T(2,:))'); % for 3D
    Volume  = norm(cross(T(1,:),T(2,:))); % for 2D
    Dens    = 2*pi/Volume/nks; 
    % Since /sigma_xy = (e^2/h)*(1/2*pi)*k-space-integral, dkx*dky = (2*pi)^2 / volume
    % B.Anderi Bernevig, Taylor L. Hughes. Topological Insulators and Topological Superconductors, eq.(3.1)
    %Ainv     = 1/norm(cross(BR(1,:),BR(2,:)));
    %Dens     = 2*pi*Ainv/N1/N2;    
    COF      = Dens*cquantum;             % integral measure of the BZ with the unit (e^2/h)
    mu       = linspace(Emu(1),Emu(2),nE);    % chemical potnetial
    Berry_ss = zeros(size(mu,2),9);           % Berry phase of each direction as function of energy
%     LM_ss    = zeros(size(mu,2),3); 

    for imu = 1:length(mu);
        u  = mu(imu);
        fn = (exp((Ek(:,:)-(Ef+u))/(kB*Temp)) + 1).^-1;     % Fermi distribution   
        for ii=1:9
%             berry_tmp{ii} = Berry(:,:,:,ii);
%             Berry_new = Berry(:,:,iband,:);
            Berry_ss(imu,ii) = COF*sum(sum(fn.*Berry(:,:,ii))); %sum over the orbitals, then, integral over k-space
%             Berry_ss(imu,ii) = COF*sum(sum(Berry(:,7,ii)));
%             LM_ss(imu,ii)    = COF*sum(sum(sum(fn.*LM(:,:,:,ii))));
        end
    end
else
    berry_k  = zeros(nks,2);
    fn = (exp((Ek(:,:)-Ef)/(kB*Temp)) + 1).^-1;
    for ii=1:9
        btmp(:,:,ii) = fn.*Berry(:,:,ii)*abc(1)*abc(2);
        berry_k    = sum(btmp,2);
    end
    figure('Name','Band','NumberTitle','off');
    plot(1:nks,Ek(:,:)-Ef);
    hold on
    line('XData',[1 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
    hold off
    
    figure('Name','SBerry','NumberTitle','off');    
    for is=1:9
        plot(1:nks,berry_k(:,1,is));
        hold on
    end
    line('XData',[1 nks], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
    return
end

%% --- Plot Figure --- %%%
%% --- Conductivity --- %
figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
       'papertype','usletter','numbertitle','off',...
       'PaperPositionMode','manual','paperorientation','portrait',...
       'color','w');
clf
for ii=1:9
    plot(mu,Berry_ss(:,ii),'Linewidth',1.0);
    hold on
end
line('XData',[mu(1) mu(end)], 'YData', [0 0], 'LineStyle', '--', ...
    'LineWidth', 0.5, 'Color','k');
% ls = legend('\bf{$\Omega_{x}$}','\bf{$\Omega_{y}$}','\bf{$\Omega_{z}$}');
% set(ls,'interpreter','LaTex','FontSize',20);
set(gca,'Ticklength',[0.05 0]);
set(gca,'Fontsize',15);
xlabel('\bf{$\mu$ (eV)}','interpreter','LaTex');
ylabel('\bf{Berry}','interpreter','LaTex');
ax = gca;
ax.TickLabelInterpreter='latex';
ax.FontSize = 20;

%--- LM ---%
% figure('position',[150 0 850 660],'paperposition',[0.25 0.25 8 10.5],...
%        'papertype','usletter','numbertitle','off',...
%        'PaperPositionMode','manual','paperorientation','portrait',...
%        'color','w');
% clf
% for ii=1:3
%     plot(mu,LM_ss(:,ii),'Linewidth',1.0);
%     hold on
% end
% line('XData',[mu(1) mu(end)], 'YData', [0 0], 'LineStyle', '--', ...
%     'LineWidth', 0.5, 'Color','k');
% ls = legend('\bf{LM$_{x}$}','\bf{LM$_{y}$}','\bf{LM$_{z}$}');
% set(ls,'interpreter','LaTex','FontSize',20);
% set(gca,'Ticklength',[0.05 0]);
% set(gca,'Fontsize',15);
% xlabel('\bf{$\mu$ (eV)}','interpreter','LaTex');
% ylabel('\bf{Berry}','interpreter','LaTex');
% ax = gca;
% ax.TickLabelInterpreter='latex';
% ax.FontSize = 20;
