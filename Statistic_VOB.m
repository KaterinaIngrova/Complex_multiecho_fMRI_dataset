%% STATISTIC FOR SINGLE SUBJECT ANALYSIS - for task VOB
function Statistic_VOB

    clear all;
    clear classes

    %Prepare and add path 
    addpath /cluster/projekty/_tools/SPM12
    rehash path;        
    spm_jobman('initcfg');
    cwd='/cluster/projekty/ELICIT/derivatives/spm_processed_dataset'; %Path to data folder

    pathdelim = '/'; %for UNIX this is '/' e.g. '/usr/home', for windows it is '\' e.g. 'c:\'

    %Define sessions - name of task, number of runs, number of echoes, indicator for phase data (0 = only mag, 1 = mag+phase)
    sessions={    
        'VOB_acq-tr1800', 1 , 3, 0
        'VOB_acq-tr800', 1 , 3, 0
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

    nsub = size(subjects,1);
    nses = size(sessions,1);
    nruns = 1;

    regr_names = {'x trans', 'y trans', 'z trans', 'x rot', 'y rot', 'z rot','dx trans', 'dy trans', 'dz trans', 'dx rot', 'dy rot', 'dz rot', 'x trans2', 'y trans2', 'z trans2', 'x rot2', 'y rot2', 'z rot2','dx trans2', 'dy trans2', 'dz trans2', 'dx rot2', 'dy rot2', 'dz rot2'}; %Name for movement regressors
    ModelDirectoryName = 'Model_24MR_ARFAST';

    spm_defaults;
    global defaults;

    %For loop across subjects
    for sub=1:nsub  
        sub_dir = fullfile(cwd,subjects{sub,1}); %Path to subject
    
        %For loop across sessions
        for ses=1:nses      
            fcnidir=fullfile(sub_dir,'func'); %Path to functional data

            %Getting TR, TE and other specific data from header
            NoOfEchoes = sessions{ses,3};
            PHflag = sessions{ses,4}; %Flag for phase images, if 1 both magnitude and phase will be processed
            direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func'];   
            for ii = 1:NoOfEchoes
                qj   = spm_select('FPList',direc{ses},['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-' num2str(ii) '_part-mag_bold\.json$']);        
                jsondata=spm_jsonread(qj);
                TEs(ii)=jsondata.EchoTime * 1000; %convert to miliseconds
                TR = jsondata.RepetitionTime;
                SliceTimes = jsondata.SliceTiming;
            end

 
            %Get preprocessed and movement files
            P   = [];
            q   = spm_select('FPList',fcnidir,['^swtse_r.*' subjects{sub,1} '.*' sessions{ses,1} '.*_part-mag_bold.nii$']);
            P   = strvcat(P,q);
            tmpv=spm_vol(P);
            nscans = size(tmpv,1);
            for i=1:nscans
                selected_scans{i,1} = [P ',' num2str(i)];
            end
    
    
            mrfile = spm_select('FPList',fcnidir,['^rp.*' subjects{sub,2} '.*' sessions{ses,1} '.*txt$']);
            xmrf= spm_load(mrfile(1,:));
            mrf = xmrf(1:nscans,:);
    
            %Load timing
            timing_file = fullfile(sub_dir, 'func', [subjects{sub,2} '_task-' sessions{ses,1} '_events.tsv']);
            TT = readtable(timing_file, 'FileType', 'text', 'Delimiter', '\t');

            stim_onset = TT.onset;
            stim_duration = TT.duration;
            EInfo.TrialType = TT.trial_type;
   
            %Find stimuli types
            EventType=zeros(1,length(EInfo.TrialType));
            for i=1:length(EInfo.TrialType)
                switch string(EInfo.TrialType(i))
                    case {'frequent'}
                        EventType(i) = 1;
                    case {'target'}
                        EventType(i) = 2;
                    case {'distractor'}
                        EventType(i) = 3;
                end
            end
   
  
            mdir = fullfile(sub_dir,[ModelDirectoryName '_' sessions{ses,1}]);
            if ~(exist(mdir))
              mkdir(mdir);
            end
            cd(mdir);
            model_dir{1}={mdir};
    
    
            timingstruct  = struct(  'units',       'secs',...
                                     'RT',          TR,...
                                     'fmri_t',      16,...
                                     'fmri_t0',     8);
    
            modstruct     = struct(  'name',        {},...
                                     'param',       {},...
                                     'poly',        {});                         
 
            %Conditions                         
            s.condstruct(1) = struct(  'name',             'Frequent',...
                                            'onset',       stim_onset(EventType==1),...
                                            'duration',    stim_duration(EventType==1),...
                                            'tmod',        {0},...
                                            'mod',         modstruct);

            s.condstruct(2) = struct(  'name',             'Target',...
                                            'onset',       stim_onset(EventType==2),...
                                            'duration',    stim_duration(EventType==2),...
                                            'tmod',        {0},...
                                            'mod',         modstruct);
                                    
            s.condstruct(3) = struct(  'name',             'Distractor',...
                                            'onset',       stim_onset(EventType==3),...
                                            'duration',    stim_duration(EventType==3),...
                                            'tmod',        {0},...
                                            'mod',         modstruct);                                    
                                                                       
       
            %Movement regressors
            moves = mrf;
            dmoves = [0 0 0 0 0 0; diff(moves)];
            fmoves =[moves dmoves moves.^2 dmoves.^2];                           
            for i = 1:24
                regressstruct(1,i) = struct(  'name',        regr_names{i},...
                                              'val',         {fmoves(1:nscans,i)});     
            end
           

                
            for ss = 1:1   
                sessstruct(ss) = struct( 'scans',       {selected_scans},...
                                     'cond',        s(ss).condstruct,...
                                     'multi',       {{''}},...
                                     'regress',     regressstruct(ss,:),...
                                     'multireg',    {{''}},...
                                     'hpf',         {160});
            end
                         
            factstruct    = struct(  'name',        {},...
                                     'levels',      {});
                         
            basesstruct.hrf.derivs = [0 0];
    
            fmrispecstruct= struct(  'dir',         model_dir,...
                                     'timing',      timingstruct,...
                                     'sess',        sessstruct,...
                                     'fact',        factstruct,...
                                     'bases',       basesstruct,...
                                     'volt',        1,...
                                     'global',      'None',...
                                     'mask',        {{''}},...
                                     'cvi',         'FAST');
                         
                            
            spm_unlink(fullfile(mdir, 'SPM.mat')); 
    
            matlabbatch{1}.spm.stats.fmri_spec = fmrispecstruct;
            output_list = spm_jobman('run',matlabbatch);
            clear matlabbatch;

    
            matfile_path{1} = {[mdir pathdelim 'SPM.mat']};

            methodstruct.Classical = 1;                         

            fmrieststruct = struct(  'spmmat',      matfile_path,...
                                     'method',      methodstruct);

            spm_unlink(fullfile(mdir, 'mask.img')); 

            matlabbatch{1}.spm.stats.fmri_est = fmrieststruct;
            spm_jobman('run',matlabbatch);
            clear matlabbatch;
            
            %Contrasts
            nmr = 24; % number of movement regressors (or other nuisance variables/covariates placed behind condition regressors)
  
            consessstruct{1}.tcon.name='Frequent';
            consessstruct{1}.tcon.convec=[1 0 0 zeros(1,nmr) 0];

            consessstruct{2}.tcon.name='Target';
            consessstruct{2}.tcon.convec=[0 1 0 zeros(1,nmr) 0];

            consessstruct{3}.tcon.name='Distractor';
            consessstruct{3}.tcon.convec=[0 0 1 zeros(1,nmr) 0];

            consessstruct{4}.tcon.name='Target - Frequent';
            consessstruct{4}.tcon.convec=[-1 1 0 zeros(1,nmr) 0];    

            consessstruct{5}.tcon.name='Target - Distractor';
            consessstruct{5}.tcon.convec=[0 1 -1 zeros(1,nmr) 0];

            consessstruct{6}.tcon.name='Distractor - Frequent';
            consessstruct{6}.tcon.convec=[-1 0 1 zeros(1,nmr) 0];     

            consessstruct{7}.tcon.name='ALL';
            consessstruct{7}.tcon.convec=[1 1 1 zeros(1,nmr) 0];


                              
            construct = struct( 'spmmat',   matfile_path,...
                                'consess',  {consessstruct});

            matlabbatch{1}.spm.stats.con = construct; 
            spm_jobman('run',matlabbatch);
            clear matlabbatch;               
            clear P selected_scans q moves mrf;
   
        end 
    end 
end