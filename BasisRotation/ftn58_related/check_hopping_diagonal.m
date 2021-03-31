function check=check_hopping_diagonal(ftn58,hop_tol)
% check that (t^mm_(\vector n))* = t^mm_(-\vector n) 
% hop_tol is the tolerance for the maximum difference of these hopping parameters

    %addpath('/home/lengentyh/Matlab/ftn58_related/prune&select_ftn58')
    %addpath('/home/lengentyh/Matlab/bandplot_eigband_buildkpts')
    %load /home/lengentyh/Matlab/kagome_debug/finerplt.mat kpt
    addpath(genpath('./ftn58_related'))
    load finerplt.mat kpt
    fprintf('start checking ... \n') 

    check = struct;    
    kp = size(kpt,1);
    kpoints = kpt(:,3:5);

    % for all orbitals
    for orb=1:ftn58(1,1)
        one_orbital = orb;

        % add column in ftn58 to classify translation
        trsl=unique(ftn58(2:end,5:7),'rows');
        [logn,ordn] = ismember(ftn58(2:end,5:7),trsl,'rows');
        tn58(:,8) = [0; ordn];

        orb2ftn58 = selectftn58(ftn58,one_orbital);
        
        % for each translation
        for n=1:size(trsl,1)
            translation = trsl(n,:);
            if translation == [0,0,0]
                continue
            end
            translation_pair = [translation;-translation];
            cell_orb2ftn58 = selectftn58_translation(orb2ftn58,translation_pair);

            % preallocations
            norb = cell_orb2ftn58(1,1);
            nbond = cell_orb2ftn58(1,2);
            H = zeros(norb,norb,kp);
            phase = zeros(kp,nbond);
            hopping = zeros(kp,nbond);
            temp=zeros(kp,nbond);
            exponent=zeros(kp,nbond);
            % construct Hamiltonian
            for k=1:kp
                kvector       = kpoints(k,1:3)*2*pi;
                Ham	      = Heigband_debug(cell_orb2ftn58,kvector);
                H(:,:,k)      = Ham.H;
                phase(k,:)    = Ham.phase;
                hopping(k,:)  = round(Ham.hopping,ceil(log10(1/hop_tol)));
                temp(k,:)     = round(Ham.temp,ceil(log10(1/hop_tol)));
                exponent(k,:) = Ham.exponent;
            end
            % check the properties in every kpoints
            ctemp 	   = temp(:,1)+temp(:,2);
            check.phase    = isequal(phase(:,1),-phase(:,2));
            check.hopping = isequal(hopping(:,1),conj(hopping(:,2)));
	    %check.hopping  = max(abs(hopping(:,1)-conj(hopping(:,2)))) < hop_tol;
            check.exponent = isequal(exponent(:,1),conj(exponent(:,2)));
            check.temp     = isreal(ctemp(:,1));
	    check.cell_orb2ftn58=cell_orb2ftn58;

            if check.phase*check.hopping*check.exponent*check.temp ~= 1
                fprintf('bug in the orbital')
                disp(orb)
                fprintf('at the translation')
                disp(translation)
                fprintf('the thing that goes wrong: \n')
                disp(check)
		fprintf('the ftn58 of this translation pair: ')
		disp(check.cell_orb2ftn58(end-1:end,4))
                break
            end
        end
    end

    fprintf('end \n');
    rmpath(genpath('./ftn58_related'))

return
