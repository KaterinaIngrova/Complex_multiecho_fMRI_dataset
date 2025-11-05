function Mafil_dicom2preBIDS_public(varargin)
% 2025/10 TS 
% This function converts DICOM files to NIFTI and organizes converted data according to BIDS(-like) format. Tailored to MAFIL CEITEC MR data - https://mafil.ceitec.cz/
% Output to empty folder only, repeated run (update/extend existing dataset) not supported.
% For large datasets it is recommended to use SSD drive paths.
% 
% BIDS information - https://bids-specification.readthedocs.io/
% BIDS validator - https://bids-standard.github.io/bids-validator/
%
% Prerequisities:  
% 1) dcm2niix (dcm2niix executable or binary has to be in same directory with this M-file) - https://github.com/rordenlab/dcm2niix/
% 2) SPM12 toolbox in Matlab path, tested with v7771 - https://www.fil.ion.ucl.ac.uk/spm/
% 3) current version limited to Windows OS only
%
% Usage:    mafil_dicom2preBIDS                              Dialog window prompt for path selection
%           mafil_dicom2preBIDS(inpath,outpath)              Specify INPUT and OUTPUT paths for batch processing

savereport=0; % option to save output of dcm2niix command window to _report_from_dcm2niix.txt, for DEBUG only as this file may contain PATHs/potentially SENSITIVE information

% % % % % % % % % % % % % % % % do not modify below this line % % % % % % % % % % % % % % % % % % % % % % % % 

code_ver='v1.19'; % code version
bids_ver='1.8.0'; % compliant BIDS version
SPM_ver='7771'; % compatible SPM version
dcm2niix_ver='v1.0.20241211'; % compatible dcm2niix version


% % % % % %      M A F I L 2 B I D S   T A B L E      % % % % % % % %                                                                                                                                          HINT: specific first, common later
%IN:1)SeriesDescription    2)FileName      3)ImageType                  4)PulseSeqDetails    OUT:5)folder    6)TASK-        7)ACQ-        8)REC-               9)DIR- (RUN)  10)ECHO-     11)FLIP-      12)INV-     13)MT-     14)PART-   15)suffix
niftiLUT={
% anatomical high-res:
 '.*t1_mpr.*'               '.*'            {'ORIGINAL';'NORM'}         '.*tfl$'                'anat'      ''               ''           'norm'               ''            ''              ''          ''         ''         ''        'T1w'    
 '.*t2-star.*'              '.*'            {'ORIGINAL'}                '.*gre$'                'anat'      ''               ''           ''                   ''            'DCM.EchoNumber' ''         ''         ''         ''        'MEGRE' 
 'R2Star_Images'            '.*'            {'DERIVED';'R2';'STAR MAP'} '.*gre$'                'anat'      ''               ''           ''                   ''            ''              ''          ''         ''         ''        'R2starmap'

% fieldmaps:
 '.*'                       '.*_e2_ph$'     {'ORIGINAL';'PHASE'}        '.*gre_field_mapping$'  'fmap'      ''               ''           ''                   ''            ''              ''          ''         ''         ''        'phasediff'  
 '.*gre_.*reference.*'      '.*_e1$'        {'ORIGINAL'}                '.*gre_field_mapping$'  'fmap'      ''               ''           ''                   ''            ''              ''          ''         ''         ''        'magnitude1'
 '.*gre_.*reference.*'      '.*_e2$'        {'ORIGINAL'}                '.*gre_field_mapping$'  'fmap'      ''               ''           ''                   ''            ''              ''          ''         ''         ''        'magnitude2'  

 % func/bold/MB+ME:
 'bold_VOB_tr800.*_SBRef$'  '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VOB'            'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_VOB_tr800.*_SBRef$'  '.*_e[1-9]$'    {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VOB'            'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  
 'bold_VOB_tr1800.*_SBRef$' '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VOB'            'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_VOB_tr1800.*_SBRef$' '.*_e[1-9]$'    {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VOB'            'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  
 'bold_VisMot_tr800.*_SBRef$' '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}      '.*ep2d_bold$'          'func'      'VisMot'         'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_VisMot_tr800.*_SBRef$' '.*_e[1-9]$'    {'ORIGINAL'}              '.*ep2d_bold$'          'func'      'VisMot'         'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  
 'bold_VisMot_tr1800.*_SBRef$' '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}     '.*ep2d_bold$'          'func'      'VisMot'         'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_VisMot_tr1800.*_SBRef$' '.*_e[1-9]$'    {'ORIGINAL'}             '.*ep2d_bold$'          'func'      'VisMot'         'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  
 'bold_REST_tr800.*_SBRef$'  '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}       '.*ep2d_bold$'          'func'      'rest'           'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_REST_tr800.*_SBRef$'  '.*_e[1-9]$'    {'ORIGINAL'}               '.*ep2d_bold$'          'func'      'rest'           'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  
 'bold_REST_tr1800.*_SBRef$' '.*_e[1-9]_ph$' {'ORIGINAL';'PHASE'}       '.*ep2d_bold$'          'func'      'rest'           'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'sbref'
 'bold_REST_tr1800.*_SBRef$' '.*_e[1-9]$'    {'ORIGINAL'}               '.*ep2d_bold$'          'func'      'rest'           'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'sbref'  

 'bold_VOB_tr800$'      '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VOB'            'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_VOB_tr800$'      '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VOB'            'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold'  
 'bold_VOB_tr1800$'     '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VOB'            'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_VOB_tr1800$'     '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VOB'            'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold'  
 'bold_VisMot_tr800$'   '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VisMot'         'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_VisMot_tr800$'   '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VisMot'         'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold'  
 'bold_VisMot_tr1800$'  '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'VisMot'         'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_VisMot_tr1800$'  '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'VisMot'         'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold'  
 'bold_REST_tr800$'     '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'rest'           'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_REST_tr800$'     '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'rest'           'tr800'      ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold'  
 'bold_REST_tr1800$'    '.*_e[1-9]_ph$'     {'ORIGINAL';'PHASE'}        '.*ep2d_bold$'          'func'      'rest'           'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'phase'   'bold'
 'bold_REST_tr1800$'    '.*_e[1-9]$'        {'ORIGINAL'}                '.*ep2d_bold$'          'func'      'rest'           'tr1800'     ''                   ''            'DCM.EchoNumber' ''         ''         ''         'mag'     'bold' 

};

% % % % % % initialization of paths and input parameters % % % % % %
if nargin>2 || nargin==1
    error('wrong number of input arguments')
end
programdir=fileparts(mfilename('fullpath')); % determine full path to dcm2niix binary/executable
if ~exist(fullfile(programdir,'dcm2niix.exe'),'file'); error('dcm2niix.exe not found, please check prerequisities!'); end
runcmd=fullfile(programdir,'dcm2niix.exe');
try
    [~,spm_ver_installed]=spm('Ver','',1); 
%     if ~strcmp(spm_ver_installed,SPM_ver) % more stringent condition check
    if str2double(spm_ver_installed) < str2double(SPM_ver) % less stringent condition check
        error(['SPM12 toolbox installed version (',spm_ver_installed,') is too old, please update'])
    end
catch
    error('SPM12 toolbox not detected in MATLAB PATH')
end


if  nargin==2
    indir=varargin{1};
    if ~exist(indir,'dir'); error('input error (path doesn''t exist)'); end % check if INPUT dir exists in batch mode 
    outdir=varargin{2};
else
    indir = uigetdir(defpth,'Select an input folder (DICOM files)'); % specify both INPUT and OUTPUT paths using dialog windows
    if indir==0, error('input directory not selected'); end    
    outdir= uigetdir(defpth,'Select an output folder (NIFTI preBIDS files)');
    if outdir==0, error('output directory not selected'); end
end

if numel(dir(outdir))>2; error('output directory is not empty!'); end % check if output dir is empty
tempdir=fullfile(outdir,'unsorted');
if ~exist(tempdir,'dir'); mkdir(tempdir); end % always true, in batch mode creates outdir if needed
tic;

% % % % % % conversion from DICOM to NIFTI first, not anomynous % % % % % % 
disp(' 0/5 plese wait, DICOM conversion (not-anonymized!) in progress ...')
if numel(dir(fullfile(indir,['**',filesep,'*.*'])))>4e4; fprintf(1,'\b detected huge number of files, you should take a coffee (or two)\n'); end
    [status,result] = system([runcmd,' -g i -l o -f sub-%g',filesep,'sub-%g_std%x_%3s_%p -ba n -o ', tempdir, ' ', indir]);
if status~=0
    fprintf(2,'conversion did not finish correctly, log saved as "_report_from_dcm2niix.txt" for debugging, some hints from log:\n\n'); 
    errpos=regexpi(result,'((?<=\n)Warning: File not large enough to store image data\w*[^\n]*)|((?<=\n)Warning: Missing images?\w*[^\n]*)|((?<=\n)Error: \w*[^\n]*)','match','all');
    if ~isempty(errpos)    
        fprintf(2, '%s\n', errpos{:})
        fprintf(2, '\n... fix error(s) or remove the above mentioned data and run conversion again!\n')
    end
    fid = fopen(fullfile(outdir,'_report_from_dcm2niix.txt'),'wt');fprintf(fid,'%s', result);fclose(fid);
    return;
end
if savereport
    fid = fopen(fullfile(outdir,'_report_from_dcm2niix.txt'),'wt');fprintf(fid,'%s', result);fclose(fid);
end
disp(' 1/5 conversion done, metadata extraction and anonymization ...')

% % % % % %  anonymization, collect metadata to metaraw variable % % % % % % 
todofiles=dir(fullfile(tempdir,['**',filesep,'*.*']));
todofiles([todofiles.isdir]) = []; % skip folders
metaraw=[]; % collect all metadata info from JSON files 
fldsanon={'PatientName','PatientID','PatientAge','PatientWeight','PatientSize','PatientBirthDate','SeriesInstanceUID','StudyInstanceUID','AcquisitionDateTime','StudyID'...
    'ImageComments','InstitutionalDepartmentName','InstitutionAddress','DeviceSerialNumber','AcquisitionTime','StationName','StudyDescription','ConsistencyInfo'}; % fields to remove for anonymization (remain: ReferringPhysicianName, AccessionNumber)

h = waitbar(0,'Processing');
for ind=1:size(todofiles,1) % go through json files, anonymize specified tags, derive age from dates and update json, store all metadata incl. paths to metaraw
    waitbar(ind/size(todofiles,1))
    fn=fullfile(todofiles(ind).folder,todofiles(ind).name);
    if any(regexpi(todofiles(ind).name,'.*.json$')) % cycle through json
        aa=spm_jsonread(fn);
        try 
            if isdatetime(datetime(aa.PatientBirthDate,'InputFormat','yyyy-MM-dd')) && isdatetime(datetime(aa.AcquisitionDateTime(1:10),'InputFormat','yyyy-MM-dd'))  % determine age when measured, save to json
                derivage=round(years(datetime(aa.AcquisitionDateTime(1:10),'InputFormat','yyyy-MM-dd')-datetime(aa.PatientBirthDate,'InputFormat','yyyy-MM-dd')));
            end
        catch
            derivage=0;
        end
        for ine=1:size(fldsanon,2) % anonymization by removing tags from json
            if isfield(aa,fldsanon{ine})
                aa=rmfield(aa,fldsanon{ine});
            end
        end
        spm_jsonwrite(fn,aa,struct('indent',' ')) % overwrite source
        aa=setfield(aa,'PatientAge',derivage); % additional metadata, do not save back to json (sex, age, path and filename)
        aa=setfield(aa,'filename',todofiles(ind).name); % store path and filename for later use
        aa=setfield(aa,'pathdir',todofiles(ind).folder);
        metaraw{end+1}=aa;
        clear aa fn
    end
end
close(h)
disp(' 2/5 anonymization done, BIDS sorting and renaming ...')
% pause

% % % % % % rename files matching to niftiLUT, add _RUN-## to avoid overwrite for repeated series, save additional json % % % % % % 
accnums=[];
metafin=[];  % all metadata (from json) in format: metafin{subject}{series}
h = waitbar(0,'Processing');
for ind=1:size(metaraw,2)
    waitbar(ind/size(metaraw,2))
    if ~isfield(metaraw{ind},'AccessionNumber') % if accnum doesnt exist, skip file (keep in unsorted)
            error(' AccessionNumber missing for %s,  MAFIL dataset? Try using override-id or override-name.\n',fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename));
    end
    m = find(strcmp(metaraw{ind}.AccessionNumber, accnums)); % gather subjects by accession_number
    if isempty(m)
        m = length(accnums)+1;  % add new entry if not found 
        accnums{m} = metaraw{ind}.AccessionNumber;
        metafin{m}=[];
    end
%     sn=metaraw{ind}.SeriesNumber;
    if isempty(metafin{m})
        csn=1;
    else
        csn=size(metafin{m},2)+1; % count new series number based on order (to properly list multiecho, mag+phase, etc. in one SN)
    end
    metafin{m}{csn}=metaraw{ind}; % sort metafin using accession_number (multiecho series are no-longer overwritten with last echo)
    
    for inff=1:size(niftiLUT,1) % try to match LUT based on both protocolname + filename (may contain echo numbers, phase images, etc.) + imagetype
        if ~isfield(metaraw{ind},'SeriesDescription') || ~isfield(metaraw{ind},'PulseSequenceDetails'); break; end
        if ~isfield(metaraw{ind},'ImageTypeText'); metaraw{ind}.ImageTypeText=metaraw{ind}.ImageType; end % FIX for old PRISMA data
        if any(regexpi(metaraw{ind}.SeriesDescription,niftiLUT{inff,1})) && any(regexpi(metaraw{ind}.filename(1:end-5),niftiLUT{inff,2})) && isequal(sum(ismember(unique([metaraw{ind}.ImageTypeText;metaraw{ind}.ImageType]),niftiLUT{inff,3})),size(niftiLUT{inff,3},1)) && any(regexpi(metaraw{ind}.PulseSequenceDetails,niftiLUT{inff,4})) && ~(metaraw{ind}.SeriesNumber > 1000 && metaraw{ind}.AcquisitionNumber > 1)   % ismember or contains
            
             % % % generate proper name
             taskname=tryextractnamestring(metaraw{ind},niftiLUT{inff,6},'_task-');
             acqname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,7},'_acq-');
             recname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,8},'_rec-');
             dirname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,9},'_dir-');
             echoname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,10},'_echo-');
             flipname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,11},'_flip-');
             invname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,12},'_inv-');
             mtname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,13},'_mt-');
             partname=tryextractnamenumeric(metaraw{ind},niftiLUT{inff,14},'_part-');
             suffixname=['_' niftiLUT{inff,15}]; % suffix is mandatory
             
             % % % special handling for FUNC (BOLD)
             if strcmp(niftiLUT{inff,5},'func') && strcmp(niftiLUT{inff,15},'bold') % in case of func, create .json file with task name, omit for SBRef and different RUNs (assume repeated sequence)
                nme=strrep(taskname,'_task-','');
                descr = struct('TaskName',nme);
                spm_jsonwrite(fullfile(outdir,['task-',nme,'_',niftiLUT{inff,15},'.json']),descr, struct('indent',' '));
             end
             
             % % % special handling for fieldmaps
             if strcmp(niftiLUT{inff,15},'phasediff') % in case of fmap phasediff, add Echotime1 and EchoTime2 to json before renaming
                 if strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-1}.SeriesDescription) && strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-2}.SeriesDescription)
                     tempjson=spm_jsonread(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename));
                     tempjson=setfield(tempjson,'EchoTime1',metaraw{ind-2}.EchoTime);
                     tempjson=setfield(tempjson,'EchoTime2',metaraw{ind-1}.EchoTime);
                     spm_jsonwrite(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename),tempjson, struct('indent',' '));
                 elseif strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-2}.SeriesDescription) && strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-3}.SeriesDescription) % CIMA.X hotfix for ND_fieldmaps
                     tempjson=spm_jsonread(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename));
                     tempjson=setfield(tempjson,'EchoTime1',metaraw{ind-3}.EchoTime);
                     tempjson=setfield(tempjson,'EchoTime2',metaraw{ind-2}.EchoTime);
                     spm_jsonwrite(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename),tempjson, struct('indent',' '));
                 elseif strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-3}.SeriesDescription) && strcmp(metaraw{ind}.SeriesDescription,metaraw{ind-4}.SeriesDescription) % CIMA.X hotfix for ND_fieldmaps
                     tempjson=spm_jsonread(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename));
                     tempjson=setfield(tempjson,'EchoTime1',metaraw{ind-4}.EchoTime);
                     tempjson=setfield(tempjson,'EchoTime2',metaraw{ind-3}.EchoTime);
                     spm_jsonwrite(fullfile(metaraw{ind}.pathdir,metaraw{ind}.filename),tempjson, struct('indent',' '));
                 else
                     fprintf(2,'warn: bogus filedmap detected, failed to get EchoTimes properly for %s_%d_%s (file kept in unsorted folder)\n',metaraw{ind}.AccessionNumber,metaraw{ind}.SeriesNumber,metaraw{ind}.SeriesDescription);
                     break % skips this file - it will remain in unsorted folder
                 end
             end
             
             % % % create folder if needed, compose name and check if exists, add _run-##, move from unsorted folder
             if ~exist(fullfile(outdir,['sub-',accnums{m}],niftiLUT{inff,5}),'dir'); mkdir(fullfile(outdir,['sub-',accnums{m}],niftiLUT{inff,5})); end
                newname=fullfile(outdir,['sub-',accnums{m}],niftiLUT{inff,5},['sub-',accnums{m},taskname,acqname,recname,dirname,echoname,flipname,invname,mtname,partname,suffixname]);
                fileorder=1;
                while exist([newname,'.json'],'file')
                    fileorder=fileorder+1; % rename to _run-##, begin with 2 and increase
                    newname=fullfile(outdir,['sub-',accnums{m}],niftiLUT{inff,5},['sub-',accnums{m},taskname,acqname,recname,dirname,'_run-',sprintf('%02.f',fileorder),echoname,flipname,invname,mtname,partname,suffixname]);
                end
                for indd={'.json','.nii','.nii.gz','.bvec','.bval'} % try all possible file extensions (generated during conversion)
                    try
                        movefile(fullfile(metaraw{ind}.pathdir,[regexprep(metaraw{ind}.filename,'.json','') indd{:}]), [newname indd{:}])
                        metafin{m}{csn}.bidsname=strrep(newname,outdir,'');   % temporary store new name for "list of sequences"
                    end
                end
             clear newname fileorder descr tempjson suffixname nme taskname acqname recname dirname echoname flipname invname mtname partname suffixname
             break 
         end % files not matching to niftiLUT will remain in unsorted
    end
end
close(h)
disp(' 3/5 BIDS sorting done, save additional BIDS files ...')

 % % % % % % dataset_description.json % % % % % % %
descr = struct('Name','Dataset measured/produced at MAFIL core facility, CEITEC, Masaryk University','BIDSVersion',bids_ver,'DatasetType','raw','GeneratedBy',{{struct('Name','mafil_dicom2preBIDS','Version',code_ver)}});
spm_jsonwrite(fullfile(outdir,'dataset_description.json'),descr, struct('indent',' '));

% % % % % % % participants.tsv - fill in using metafin json information % % % % % % %
participants=[];
reportseq=[];
cntr=1;
for ing=1:size(metafin,2)
    participants.participant_id{ing}=(['sub-',metafin{ing}{1}.AccessionNumber]);
    try
        participants.sex{ing}=metafin{ing}{1}.PatientSex;
    catch
        participants.sex{ing}='o';
    end
    participants.age(ing)=metafin{ing}{1}.PatientAge;
    % participants.project{ing}=metafin{ing}{1}.ReferringPhysicianName;
    metafin{ing}(cellfun(@isempty, metafin{ing})) = []; % remove all empty cell
    for inh=1:numel(metafin{ing})
        reportseq.participant_id{cntr}=(['sub-',metafin{ing}{1}.AccessionNumber]);
        reportseq.series_number{cntr}=num2str(metafin{ing}{inh}.SeriesNumber);
        if ~isfield(metafin{ing}{inh},'SeriesDescription')
            metafin{ing}{inh}.SeriesDescription='';
            fprintf(2,'WARNING, %s missing SeriesDescription tag\n',metafin{ing}{inh}.filename);
        end
        reportseq.series_description_and_imagetype{cntr}=([metafin{ing}{inh}.SeriesDescription,regexp(metafin{ing}{inh}.filename,'(_c[0-9]*|_i[0-9]*|_e[1-9]*_ph|_e[1-9]*|_ph|_imaginary|_real|_phMag|_MoCo)(?=.json$)','match','once'),' (',metafin{ing}{inh}.ImageType{1},')']); % based on https://github.com/rordenlab/dcm2niix/blob/master/FILENAMING.md    % add file suffix from conversion (_ph, _e1) to SeriesDescription
        if isfield(metafin{ing}{inh},'bidsname')
            reportseq.BIDS_rename{cntr}=([' => ',metafin{ing}{inh}.bidsname]);
        else
            if regexpi(metafin{ing}{inh}.filename,('.*localizer.*|.*AAHead_Scout.*|.*FastView.*|.*loc_.*'),'once')  % delete unuseful files such as localizers, keep tied with cleanup section below!
                reportseq.BIDS_rename{cntr}=('x');
            elseif any(ismember(metafin{ing}{inh}.ImageType,'DERIVED')) && any(ismember(metafin{ing}{inh}.ImageType,{'MPR','CSA MPR','CSAPARALLEL'})) % delete derived images based on ImageType filter below, keep tied with cleanup section below!
                reportseq.BIDS_rename{cntr}=('x');
            else
                reportseq.BIDS_rename{cntr}=('-');
            end
        end
        cntr=cntr+1;
    end
end
inm=fieldnames(reportseq);
for inn=1:numel(inm); reportseq.(inm{inn}){cntr}=('...'); end % add dummy text to bottom line
reportseq.series_number{cntr}=(' HINT: localizers and reconstructions are deleted (x), other derived + not-recognized are kept in unsorted folder (-)'); % add hint to bottom line
if ~isempty(participants);spm_save(fullfile(outdir,'participants.tsv'),participants);end  
if ~isempty(reportseq);spm_save(fullfile(outdir,'_list_of_original_sequences.tsv'),reportseq);end

% % % % % % % participants.json % % % % % % %
jsondescr = struct('sex',struct('Description','Sex of the participant (DICOM derived)','levels',struct('M','male','F','female','o','unidentified')),'age',struct('Description','Age of the participant (DICOM derived)','Units','years'));  % from DICOM info
spm_jsonwrite(fullfile(outdir,'participants.json'),jsondescr, struct('indent',' '));

% % % % % % % default README and .bidsignore files % % % % % % %
fid = fopen(fullfile(outdir,'.bidsignore'),'wt');fprintf(fid,'unsorted/\n*.mat\n*.url\n_list_of_original_sequences.tsv\n_report_from_dcm2niix.txt');fclose(fid);
fid = fopen(fullfile(outdir,'README'),'wt');fprintf(fid,['RAW MR data converted to NIFTI+JSON(metadata)\nTo be filled in...\n\n']);fclose(fid);

disp(' 4/5 BIDS metadata saved, cleaning "unsorted" folder ...')

% % % % % % % CLEANUP % % % % % % %
toremove=dir(fullfile(tempdir,['**',filesep,'*.*']));
idx = cellfun(@isempty,regexpi({toremove.name},{'.*localizer.*|.*AAHead_Scout.*|.*FastView.*|.*loc_.*'},'once'));   % delete unuseful files such as localizers, keep tied with reportseq section above!
toremove(idx)=[];
for k=1:numel(toremove)
    delete(fullfile(toremove(k).folder,toremove(k).name))
end

toremovederived=dir(fullfile(tempdir,['**',filesep,'*.json']));     % delete derived images based on ImageType filter below, keep tied with reportseq section above!
for m=1:numel(toremovederived)
    aa=spm_jsonread(fullfile(toremovederived(m).folder,toremovederived(m).name));
    if any(ismember(aa.ImageType,'DERIVED')) && any(ismember(aa.ImageType,{'MPR','CSA MPR','CSAPARALLEL'}))
        for extind={'.json','.nii','.nii.gz','.bvec','.bval'} % try all possible file extensions
            if exist(regexprep(fullfile(toremovederived(m).folder,toremovederived(m).name),'.json',extind),'file')
                delete(regexprep(fullfile(toremovederived(m).folder,toremovederived(m).name),'.json',extind));
            end
        end
    end
end

toremovedir=dir(fullfile(tempdir,['**',filesep,'*.*']));
toremovedir(~[toremovedir.isdir])=[];        % select directories only
for k=1:numel(toremovedir)  
    try
        rmdir(fullfile(toremovedir(k).folder,toremovedir(k).name))  % remove empty directories only
    end
end
% disp(' 5/5 successfully done')
t=toc;
fprintf(' 5/5 successfully done, elapsed time %1.0fs (%1.1fh)\n',t,t/3600);
if exist(tempdir,'dir')  % check if there any any not matched files (not found at niftiLUT) left and notify user
    disp('Check "unsorted" folder for possibly useful data!');
    fid = fopen(fullfile(outdir,'README'),'at');fprintf(fid,'Check "unsorted" folder for possibly useful data.');fclose(fid);
end
return

function val=tryextractnamenumeric(s,fld,idname)
val=[];
if regexp(fld,'DCM.')
    strt=regexp(fld,'DCM.');
    if isfield(s,fld(strt+4:end))
        val=[idname fld(1:strt-1) num2str(s.(fld(strt+4:end)))];
    end
elseif ~isempty(fld)
    val=[idname fld];
end    

function val=tryextractnamestring(s,fld,idname)
val=[];
if regexp(fld,'DCM.')
    strt=regexp(fld,'DCM.');
    if isfield(s,fld(strt+4:end))
        val=[idname fld(1:strt-1) strrep(matlab.lang.makeValidName(s.(strrep(fld,'DCM.','')),'ReplacementStyle','delete'),'_','')];
    end
elseif ~isempty(fld)
    val=[idname fld];
end    
