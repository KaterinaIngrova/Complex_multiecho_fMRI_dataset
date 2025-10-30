%% PREPROCESS MULTI ECHO DATA

function Preprocess_MultiEcho

    clear all;
    clear classes

    %Prepare and add path 
    addpath /cluster/projekty/_tools/SPM12 
    rehash path;        
    spm_jobman('initcfg');
    cwd='/cluster/projekty/ELICIT/derivatives/spm_processed_dataset'; %Path to data´s folder

    pathdelim = '/'; %for UNIX this is '/' e.g. '/usr/home', for windows it is '\' e.g. 'c:\'

    %Define sessions - name of task, number of runs, number of echoes, indicator for phase data (0 = only mag, 1 = mag+phase)
    sessions={    
        'rest_acq-tr1800', 1 , 3, 0
        'rest_acq-tr800', 1 , 3, 0
        'VOB_acq-tr1800', 1 , 3, 0
        'VOB_acq-tr800', 1 , 3, 0
        'VisMot_acq-tr1800', 1 , 3, 0
        'VisMot_acq-tr800', 1 , 3, 0
        };


    TEs = [];

    %Define subject designation
    subjects={
        'sub-a01'
        'sub-a02' 
        'sub-a03' 
        'sub-a04'
        'sub-a05' 
        'sub-a06' 
        'sub-a07' 
        'sub-a08' 
        'sub-a09' 
        'sub-a10' 
        'sub-a11' 
        'sub-a12'
        'sub-a13'
        'sub-a14' 
        'sub-a15' 
        'sub-a16' 
        'sub-a17' 
        'sub-a18'
        'sub-a19' 
        'sub-a20' 
        'sub-a21' 
        'sub-a22' 
        'sub-a23' 
        'sub-a24' 
        'sub-a25' 
        'sub-a26' 
        'sub-a27' 
        'sub-a28' 
        'sub-a29'
        'sub-a30'
        'sub-a31'
        'sub-a32' 
        'sub-a33' 
        'sub-a34' 
        'sub-a35'
        'sub-a36' 
        'sub-a37'
        'sub-a38' 
        'sub-a39' 
        'sub-a40' 
        'sub-a41' 
        'sub-b42' 
        'sub-b43' 
        'sub-b44' 
        'sub-b45' 
        'sub-b46' 
        'sub-b47' 
        'sub-b48' 
        'sub-b49' 
        'sub-b50' 
        'sub-b51' 
        'sub-b52'
        'sub-b53' 
        'sub-b54'
        'sub-b55'
        'sub-b56' 
        'sub-b57' 
        'sub-b58' 
        'sub-b59' 
        'sub-b60'
        'sub-b61' 
        'sub-b62' 
        'sub-b63'
        'sub-b64'
        'sub-b65'
        'sub-b66'
        'sub-b67'
        'sub-b68' 
        'sub-b69'
        'sub-b70'
        'sub-b71' 
        'sub-b72'
        'sub-b73'
        'sub-b74'
        'sub-b75'
        'sub-b76'
        'sub-b77'
        'sub-b78' 
        'sub-b79' 
        'sub-b80'
        'sub-b81'
        'sub-b82' 
        'sub-b83' 
    };

    if size(subjects,2)==1
        subjects(:,2)=subjects(:,1);
    end

    nsub=size(subjects,1);
    nses=size(sessions,1);

    reslicevoxelsize = [3 3 3]; %Voxel-size after normalization
    boundingbox = [-78 -112 -50; 78 76 85]; %FOV after normalization
    SmoothFWHMmm = [5 5 5]; %FWHM of Gaussian Kernel in 3D

    spm_defaults;
    global defaults;

    %For loop across subjects
    for sub=1:nsub  
        %For loop across sessions
        for ses=1:nses
            swd=[cwd pathdelim subjects{sub,1}]; %Subject specific directory
            cd(cwd)

            %Getting TR, TE and other specific data from header
            NoOfEchoes = sessions{ses,3};
            PHflag = sessions{ses,4}; %Flag for phase images, if 1 both magnitude and phase will be processed
            direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func']; %Path to functional data
            for ii = 1:NoOfEchoes
                qj   = spm_select('FPList',direc{ses},['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-' num2str(ii) '_part-mag_bold\.json$']);        
                jsondata=spm_jsonread(qj);
                TEs(ii)=jsondata.EchoTime * 1000; %Convert to miliseconds
                TR = jsondata.RepetitionTime;
                SliceTimes = jsondata.SliceTiming;
            end


            %-----------------------------------------------------------------    
            % Calculate mean image
            %-----------------------------------------------------------------
            %Get files to realign
            selected_scans = [];
            direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func']; %Path to functional data
            q   = spm_select('FPList',direc{ses},['^tse_r' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_bold\.nii$']);        

            tmpv = spm_vol(q);
            tmpY = spm_read_vols(tmpv);
            meanY = mean(tmpY,4);

            meanvol = tmpv(1);
            fname=meanvol.fname;
            [pathstr,name,ext]=fileparts(fname);
            meanvol.fname=fullfile(pathstr,['mean' name ext]);
            meanvol.descrip=['mean image'];
            meanvol=spm_write_vol(meanvol,meanY);

            %-----------------------------------------------------------------    
            % Set up Coregister - mean BOLD and SBRef from first echo
            %-----------------------------------------------------------------
            %Get reference image
            dirref = [cwd pathdelim subjects{sub,1} pathdelim 'func']; 
            refimg{1} = spm_select('FPList',dirref,['^meantse.*' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_bold\.nii$']);

            %Get source image 
            dirsrc = [cwd pathdelim subjects{sub,1} pathdelim 'func']; 
            srcimg{1} = spm_select('FPList',dirsrc,['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_sbref\.nii$']);

            %Other images (nothing in this script)
            othimg = {''};

            COR.estimate.ref = refimg;
            COR.estimate.source = srcimg;
            COR.estimate.other = othimg;
            COR.estimate.eoptions = struct( 'cost_fun',   'nmi',...
                                            'sep',        [4 2],...
                                            'tol',        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001],...
                                            'fwhm',       [7 7]);

            matlabbatch{1}.spm.spatial.coreg = COR;                                 
            spm_jobman('run',matlabbatch);
            clear matlabbatch COR dirsrc refimg srcimg dirref;    

            %-----------------------------------------------------------------    
            % Spatial normalisation SBRef
            %-----------------------------------------------------------------
            dirsrc = [cwd pathdelim subjects{sub,1} pathdelim 'func']; 
            srcimg{1} = spm_select('FPList',dirsrc,['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_sbref\.nii$']);

            SUBstr.vol = srcimg;
            SUBstr.resample = srcimg;                                           

            EOstr = struct ('biasreg',  0.0001,...
                            'biasfwhm', 60,...
                            'tpm',      {[spm('dir') '\tpm\TPM.nii']},...
                            'affreg',   'mni',...
                            'reg',      [0 0.001 0.5 0.05 0.2],...
                            'fwhm',     0,...
                            'samp',     3);

            WOstr.bb = [-78 -112 -70; 78 76 85];
            WOstr.vox = [3 3 3];
            WOstr.interp = 4;

            ESTWR = struct ('subj',     SUBstr,...
                            'eoptions', EOstr,...
                            'woptions', WOstr);

            matlabbatch{1}.spm.spatial.normalise.estwrite = ESTWR;
            spm_jobman('run',matlabbatch);
            clear matlabbatch SUBstr WOstr EOstr ESTWR srcimg;

            %-----------------------------------------------------------------    
            % Spatial normalisation fMRI
            %-----------------------------------------------------------------
            dirsrc = [cwd pathdelim subjects{sub,1} pathdelim 'func']; 
            srcimg{1} = spm_select('FPList',dirsrc,['^y_' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_sbref\.nii$']);

            SUBstr.def = srcimg;
            tmpf=[];
            tmpf = strvcat(tmpf, spm_select('FPList',fullfile(swd,'func'),['^tse.*' subjects{sub,2} '.*' sessions{ses,1} '_echo-1_part-mag_bold\.nii$']));
            tmpf = strvcat(tmpf, spm_select('FPList',fullfile(swd,'func'),['^meantse_r.*' subjects{sub,2} '.*' sessions{ses,1} '_echo-1_part-mag_bold\.nii$']));
            resamp=cellstr(tmpf);

            SUBstr.resample = resamp;

            WOstr.bb = boundingbox;
            WOstr.vox = reslicevoxelsize;
            WOstr.interp = 4;

            WRITE = struct ('subj',     SUBstr,...
                            'woptions', WOstr);

            matlabbatch{1}.spm.spatial.normalise.write = WRITE;
            spm_jobman('run',matlabbatch);
            clear matlabbatch SUBstr WOstr EOstr tmpf WRITE dirsrc srcimg;


            %-----------------------------------------------------------------    
            % Set up Smooth functional images
            %-----------------------------------------------------------------
            %Get files to smooth
            P     = [];
            direc  = fullfile(swd,'func');  
            q   = spm_select('FPList',direc,['^wtse_r.*' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-1_part-mag_bold\.nii$']);
            P   = strvcat(P,q);

            for i=1:size(P,1)
                selected_scans{i,1} = [P(i,:)];
            end

            SMT.data = selected_scans;
            SMT.fwhm = SmoothFWHMmm;
            SMT.dtype = 0;
            SMT.im = 0;
            SMT.prefix = 's';

            matlabbatch{1}.spm.spatial.smooth = SMT;
            spm_jobman('run',matlabbatch);
            clear matlabbatch;

        clear P dirref direc SMT dirsrc srcimg refimg othimg COR RW q matlabbatch selected_scans RW COR NRA NRF WRITE WOstr SUBstr resamp;

        end 
    end 
end
