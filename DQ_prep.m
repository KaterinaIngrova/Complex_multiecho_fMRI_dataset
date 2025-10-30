%% COMPUTING DATA QUALITY MATRICS (tSNR map, SNS map, DVARS, FD)

function DQ_prep;

    %Prepare and add path     
    addpath '/cluster/projekty/_tools/SPM12';
    addpath '/cluster/projekty/ELICIT/derivatives/spm_processed_dataset/code'; %Folder with all codes
    rehash path; 
    spm_jobman('initcfg');
    DefStruc.StudyDir='/cluster/projekty/ELICIT/derivatives/spm_processed_dataset';  %Path to data´s folder
    cwd=DefStruc.StudyDir;

    DefStruc.pathdelim = filesep; %Platform specific delimiter in path, e.g. / or \

    %Define sessions - name of task, number of runs, number of echoes, indicator for phase data (0 = only mag, 1 = mag+phase)
    sessions={    
        'rest_acq-tr1800', 1 , 3, 0
        'rest_acq-tr800', 1 , 3, 0
        'VOB_acq-tr1800', 1 , 3, 0
        'VOB_acq-tr800', 1 , 3, 0
        'VisMot_acq-tr1800', 1 , 3, 0
        'VisMot_acq-tr800', 1 , 3, 0
        };
    
    %Define subject designation
    DefStruc.subjects={
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
  
    if size(DefStruc.subjects,2)==1
        DefStruc.subjects(:,2)=DefStruc.subjects(:,1);
    end

    nsub=size(DefStruc.subjects,1);
    nses=size(DefStruc.sessions,1);

    DefStruc.mask_dir = '/cluster/projekty/ELICIT/derivatives/spm_processed_dataset/masks'; %Directory with ROI and other masks
    DefStruc.WMmask_filename =  fullfile(DefStruc.mask_dir,'white.nii');
    DefStruc.GMmask_filename = fullfile(DefStruc.mask_dir,'grey.nii');
    DefStruc.CSFmask_filename = fullfile(DefStruc.mask_dir,'csf.nii');
    DefStruc.FuncNetworks_filename = fullfile(DefStruc.mask_dir,'aal_2_networks.mat');
    DefStruc.FuncNetMask_filename = fullfile(DefStruc.mask_dir,'aal.nii');
    DefStruc.WMmask_threshold = 0.9;
    DefStruc.GMmask_threshold = 0.8;
    DefStruc.CSFmask_threshold = 0.9;
    DefStruc.datafileprep = 'wtse_r'; %Preposition of files with preporcesed data 
    DefStruc.datafileext = 'nii';
    DefStruc.ibv_threshold = 0.5;
    DefStruc.obv_threshold = 0.1;

    %For loop across subjects
    for sub=1:nsub
        %For loop across sessions
        for ses =1:nses
            datadir = fullfile(DefStruc.StudyDir,DefStruc.subjects{sub,1},'func'); %Path to functional data
            inputdata = spm_select('FPList',datadir,['^' DefStruc.datafileprep DefStruc.subjects{sub,2} '.*' DefStruc.sessions{ses,1} '.*' DefStruc.datafileext '$']);
            
            %Getting TR, TE and other specific data from header
            for ii = 1
                qj   = spm_select('FPList',datadir,['^' DefStruc.subjects{sub,2} '.*' DefStruc.sessions{ses,1} '_.*echo-' num2str(ii) '_part-mag_bold\.json$']);        
                jsondata=spm_jsonread(qj);
                TEs(ii)=jsondata.EchoTime * 1000; %Convert to miliseconds
                TR = jsondata.RepetitionTime;
                SliceTimes = jsondata.SliceTiming;
            end

            DefStruc.TR = TR;

            DQ_calc(inputdata,DefStruc,0);

            %Calcluate FD from realignment parameters
            qmp   = spm_select('FPList',datadir,['^rp_' DefStruc.subjects{sub,2} '.*' DefStruc.sessions{ses,1} '_.*echo-2_part-mag_bold\.txt$']);  
            mp=spm_load(qmp);
            FD=FD_calc(mp);
            fdfilename=fullfile(datadir,'DQoutputs',['FDout_' DefStruc.subjects{sub,2} '_' DefStruc.sessions{ses,1} '.mat']);
            save(fdfilename,'FD');
        end 
    end 
end 


