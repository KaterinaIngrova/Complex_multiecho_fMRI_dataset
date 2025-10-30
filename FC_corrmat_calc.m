%% COMPUTATION OF CORRELATION MATRICES AND SIGNAL QUALITY METRICS 
function FC_corrmat_calc;

    %Prepare and add path       
    addpath /cluster/projekty/_tools/SPM12
    rehash path;        
    spm_jobman('initcfg');
    DefStruc.StudyDir='/cluster/projekty/ELICIT/derivatives/spm_processed_dataset/'; %Path to data´s folder
    cwd=DefStruc.StudyDir;

    DefStruc.pathdelim = filesep; %Platform specific delimiter in path, e.g. / or \

    %Define session - name of task, number of runs, number of echoes, indicator for phase data (0 = only mag, 1 = mag+phase)
    DefStruc.sessions={    
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

    nsub=size(DefStruc.subjects,1);
    nses=size(DefStruc.sessions,1);

    DefStruc.ModelDirectoryName = 'FC_CorrMatrix';
    DefStruc.mask_dir = '/cluster/projekty/ELICIT/derivatives/spm_processed_dataset/masks'; %Directory with ROI and other masks
    DefStruc.ROI_filename = 'aal.nii'; %Filename of ROI mask
    DefStruc.ROInameshort = 'AAL_tse'; %Shor name/description of ROIs used for naming of directories with results
    DefStruc.WMmask_filename = 'white.nii';
    DefStruc.GMmask_filename = 'grey.nii';
    DefStruc.CSFmask_filename = 'csf.nii';
    DefStruc.FuncNetworks_filename = 'aal_2_networks.mat';
    DefStruc.WMmask_threshold = 0.9;
    DefStruc.GMmask_threshold = 0.8;
    DefStruc.CSFmask_threshold = 0.9;
    DefStruc.datafileprep = 'wtse_r.*echo-1'; %Preposition of files with preporcesed data 
    DefStruc.datafileext = 'nii';
    DefStruc.ibv_threshold = 0.8;
    DefStruc.obv_threshold = 0.1;

    %For loop across subjects
    for sub=1:nsub  
        %For loop across sessions
        for ses =1:nses % 
            ROI_IDs = [1:116]; %(ROI_ids within one mask file used for correlation matrix)

            [RepreSigs, ROIs_info, DefStruc] = Representatives_routine(DefStruc,sub,ses,ROI_IDs); 

            Results{sub,ses}.RepreSigs = RepreSigs;
            Results{sub,ses}.ROIs_info = ROIs_info;

            [CorrMat CorrMatFull] = PCorrelationMatrix(0.2, RepreSigs, ROIs_info, DefStruc);

            Results{sub,ses}.CorrMat = CorrMat;
            Results{sub,ses}.CorrMatFull = CorrMatFull;

            %Save single subject data
            savedir = fullfile(DefStruc.StudyDir,DefStruc.subjects{sub},DefStruc.ModelDirectoryName);
            if ~(exist(savedir))
                mkdir(savedir);             
            end
            display(savedir);
            savefilename = fullfile(savedir,['FCmatrix_' DefStruc.ROInameshort '_' DefStruc.subjects{sub,1} '_' DefStruc.sessions{ses,1}]);
            save(savefilename,'RepreSigs','ROIs_info','CorrMat','CorrMatFull');
        end        
    end 

    savefilename = fullfile(DefStruc.StudyDir,['FCmatrix_' DefStruc.ROInameshort '.mat']);
    save(savefilename,'Results','DefStruc');
end 



function [RepreSigs, ROIs_info, DefStruc]=Representatives_routine(DefStruc,sub,ses,ROI_IDs) 
    %Create model directory
    mdir = fullfile(DefStruc.StudyDir,DefStruc.subjects{sub,1});
    
    %Prepare confounds for data filtering
    datadir = fullfile(DefStruc.StudyDir,DefStruc.subjects{sub,1},'func');
    
    %Movement
    mov_txt = spm_select('FPList', datadir, ['^rp_.*' DefStruc.sessions{ses,1} '.*\.txt$']);
    mov_par=load(mov_txt);
    mov_par2=mov_par.^2;
    mov_par_diff=[zeros(1,size(mov_par,2));diff(mov_par)]; %First row are zeros, differences are following
    mov_par_diff2=mov_par_diff.^2;
    mrfx=[mov_par mov_par2 mov_par_diff mov_par_diff2];
    
    %WM and CSF confounds - preparation
    %Resample all masks into the same space accordint to first data file
    wmmask = spm_select('FPList',DefStruc.mask_dir,['^' DefStruc.WMmask_filename '$']);
    gmmask = spm_select('FPList',DefStruc.mask_dir,['^' DefStruc.GMmask_filename '$']);
    csfmask = spm_select('FPList',DefStruc.mask_dir,['^' DefStruc.CSFmask_filename '$']);
    roimask = spm_select('FPList',DefStruc.mask_dir,['^' DefStruc.ROI_filename '$']);
    
    mov_txt = spm_select('FPList', datadir, ['^rp_.*' DefStruc.sessions{ses,1} '.*\.txt$']);
    datafiles = spm_select('FPList',datadir,['^swtse_rsub-.*' DefStruc.sessions{ses,1} '.*\.' DefStruc.datafileext '$']);
    disp('Found data files:')
    disp(datafiles)
    
    ROIs_info.ROImaskfile = DefStruc.ROI_filename;
    
    P = {
     [datafiles(1,:) ',1']
     wmmask
     csfmask
     gmmask
     roimask
     };

    flags.mask=0;
    flags.mean=0;
    flags.interp=0;
    flags.which=1;
    flags.wrap=[0 0 0];
    flags.prefix = 'q';
    
    spm_reslice(P,flags)
    
    [pathstr,name,ext]=fileparts(wmmask);
    wmmask=fullfile(pathstr,['q' name ext]);
    [pathstr,name,ext]=fileparts(gmmask);
    gmmask=fullfile(pathstr,['q' name ext]);
    [pathstr,name,ext]=fileparts(csfmask);
    csfmask=fullfile(pathstr,['q' name ext]);
    [pathstr,name,ext]=fileparts(roimask);
    roimask=fullfile(pathstr,['q' name ext]);
    
    %Read all necessary data files incuding masks/ROIs
    wmv=spm_vol(wmmask);
    [wmY, ~]=spm_read_vols(wmv);
    gmv=spm_vol(gmmask);
    [gmY, ~]=spm_read_vols(gmv);
    csfv=spm_vol(csfmask);
    [csfY, ~]=spm_read_vols(csfv);
    roiv=spm_vol(roimask);
    [roiY, ~]=spm_read_vols(roiv);
    datav=spm_vol(datafiles);
    [dataY, ~]=spm_read_vols(datav);
    
        
    %Find valid (inbrain) voxels and noise (background) voxel  (subject specific mask)
    nscans = size(datav,1);
    reduced_nscans = nscans; %Or smaller number as 100 for faster comupatition
    ibv_ii=zeros(size(dataY(:,:,:,1)));
    obv_ii=zeros(size(dataY(:,:,:,1)));
    for scan=1:min(nscans,reduced_nscans)
        Yi = dataY(:,:,:,scan);
        meanY = mean(Yi(~isnan(Yi)));
        ibvthr=meanY*DefStruc.ibv_threshold;
        ibv_ii = ibv_ii + ((Yi>ibvthr)&(~isnan(Yi)));
        obvthr=meanY*DefStruc.obv_threshold;
        obv_ii = obv_ii + ((Yi<obvthr)&(~isnan(Yi)));
    end
    valid_indices = (ibv_ii==min(nscans,reduced_nscans)); %mask/indices of valid voxels 
    validindvec = find(valid_indices);
    noise_indices = (obv_ii==min(nscans,reduced_nscans)); %mask/indices of noise voxels 
    noiseindvec = find(noise_indices);
    
    dataY = reshape(dataY,[],length(datav))'; %Reshape 4D data to 2D
    
    WM_indices = (wmY>=DefStruc.WMmask_threshold)&valid_indices;
    CSF_indices = (csfY>=DefStruc.CSFmask_threshold)&valid_indices;
    GM_indices = (gmY>=DefStruc.GMmask_threshold)&valid_indices;
    WMindvec=find(WM_indices);
    GMindvec=find(GM_indices);
    CSFindvec=find(CSF_indices);
    

    %Remove constant signal to avoid NaN in subsequent calculations
    temp_noise_data = dataY(:,noiseindvec);   
    ndstd = std(temp_noise_data);
    new_noiseindvec = noiseindvec(ndstd~=0);
    noiseindvec = new_noiseindvec;
    noise_data = temp_noise_data(:,ndstd~=0); %Final noise data
    
    noise_datad = spm_detrend(noise_data,1); %Remove mean and linear drift
    meannoise = mean(mean(noise_data,1));
    stdnoise = mean(std(noise_data,0,2));
    stdnoiseTS = mean(std(noise_data,0,1));    
    
    
    GloSIg = spm_global(datav);
    GloMean = mean(mean(dataY(:,validindvec)));
    
    CSF_data = dataY(:,CSFindvec);
    CSF_datad = spm_detrend(CSF_data,1); %Remove mean and linear drift
    meanCSF = mean(CSF_data,1);
    stdCSF = std(CSF_data,0,2);
    stdCSFTS = std(CSF_data,0,1);    
    WM_data = mm_detrend(dataY(:,WMindvec),1);
    meanWM = mean(WM_data,1);
    stdWM = std(WM_data,0,2);
    stdWMTS = std(WM_data,0,1);
    tSNR_WM = meanWM./stdWMTS*sqrt(length(datav)); %Calculate tSNR
    tSNR_WMm = mean(tSNR_WM);
    tSNR_WMs = std(tSNR_WM);
    WM_datad = spm_detrend(WM_data,1); %Remove mean and linear drift
    GM_data = mm_detrend(dataY(:,GMindvec),1);
    meanGM = mean(GM_data,1);
    stdGM = std(GM_data,0,2);
    stdGMTS = std(GM_data,0,1);
    tSNR_GM = meanGM./stdGMTS*sqrt(length(datav)); %Calculate tSNR
    tSNR_GMm = mean(tSNR_GM);   
    tSNR_GMs = std(tSNR_GM); 
    GM_datad = spm_detrend(GM_data,1); %Remove mean and linear drift

    SNR_GM = mean(meanGM(:)) / meannoise; 
    SNR_WM = mean(meanWM(:)) / meannoise; 
    SNR_GM2 = mean(meanGM(:)) / stdnoise; 
    SNR_WM2 = mean(meanWM(:)) / stdnoise; 
    SNR_GM3 = mean(meanGM(:)) / (stdnoise/mean(stdGM));
    SNR_WM3 = mean(meanWM(:)) / (stdnoise/mean(stdWM));
    
    %Prepare K matrix for high-pass filtering
    for ii = 1 %Just from the first echo
   		qj   = spm_select('FPList',datadir,['^sub-.*' DefStruc.sessions{ses,1} '_.*echo-' num2str(ii) '_part-mag_bold\.json$']); 
        display(qj);
        jsondata=spm_jsonread(qj);
        TEs(ii)=jsondata.EchoTime * 1000; %Convert to miliseconds
        TR = jsondata.RepetitionTime;
        SliceTimes = jsondata.SliceTiming;
    end
    
    
    K.RT = TR;
    K.row = 1:nscans;
    K.HParam = 128;
    K = spm_filter(K);
    
    %Remove low freequences from nuisance source time-series
    WM_datad = spm_filter(K, WM_datad);
    CSF_datad = spm_filter(K, CSF_datad);

    %Representative signals as first principal component
    [WM_components WMvar] = PCA_MG(WM_datad);
    WMrepre = WM_components(:,1);
    [CSF_components CSFvar] = PCA_MG(CSF_datad);
    CSFrepre = CSF_components(:,1);

    noise_dataf = spm_filter(K, noise_datad);
    CSF_dataf = spm_filter(K, CSF_datad);
    noise_dataf2 = noise_dataf;
    [noise_components noisevar] = PCA_MG(noise_dataf);
    NOISErepre = noise_components(:,1);
    
    %Prepare Xfilt matrix for GLM based filtering  
    mrf=mrfx(1:nscans,:);
    tmpX = [mrf WMrepre CSFrepre ones(nscans,1)];    
    Xfilt = spm_filter(K, tmpX);
    
    %Filter noise signal
    noise_dataf = noise_dataf-Xfilt*(pinv(Xfilt)*noise_dataf);
    CSF_dataf = CSF_dataf-Xfilt*(pinv(Xfilt)*CSF_dataf);
    [noise_componentsf noisevarf] = PCA_MG(noise_dataf);
    stdnoiseTSF = mean(std(noise_dataf,0,1));   
    stdCSFTSF = mean(std(CSF_dataf,0,1)); 
    
    SFS_WM1 = (meanWM./GloMean).*(stdWMTS./stdnoiseTSF);
    SFS_WM1m = mean(SFS_WM1);
    SFS_WM1s = std(SFS_WM1);
    SFS_WM2 = (meanWM./GloMean).*(stdWMTS./stdCSFTSF);
    SFS_WM2m = mean(SFS_WM2);
    SFS_WM2s = std(SFS_WM2);
    SFS_GM1 = (meanGM./GloMean).*(stdGMTS./stdnoiseTSF);
    SFS_GM1m = mean(SFS_GM1);
    SFS_GM1s = std(SFS_GM1);
    SFS_GM2 = (meanGM./GloMean).*(stdGMTS./stdCSFTSF);
    SFS_GM2m = mean(SFS_GM2);
    SFS_GM2s = std(SFS_GM2); 
    
    %Extract ROI signals    
    for roi = 1:length(ROI_IDs)    
        ROI_id = ROI_IDs(roi);
        ROI_indices = (roiY==ROI_id)&valid_indices; %Mask/indices of actual ROI_id
        ROIindvec = find(ROI_indices);
    
        %Extract signal only for existing valid voxels
        if ROIindvec 
            ROI_data = dataY(:,ROIindvec);
          
            %Calculate SNR and tSNR for each ROI
            ROImean = mean(ROI_data,1);
            ROIstd = std(ROI_data,0,1);
            tSNR = ROImean./ROIstd;
            ROIs_info.SNR{roi} = mean(ROImean(:))/meannoise;
            ROIs_info.SNR2{roi} = mean(ROImean(:))/stdnoise;
            ROIs_info.SNR3{roi} = mean(ROImean(:))/(stdnoise/mean(ROIstd));
            ROIs_info.tSNRm{roi} = mean(tSNR);
            ROIs_info.tSNRs{roi} = std(tSNR);
            SFS1 = (ROImean./GloMean).*(ROIstd./stdnoiseTSF);
            ROIs_info.SFS1m{roi} = mean(SFS1);
            ROIs_info.SFS1s{roi} = std(SFS1);
            SFS2 = (ROImean./GloMean).*(ROIstd./stdCSFTSF);
            ROIs_info.SFS2m{roi} = mean(SFS2);
            ROIs_info.SFS2s{roi} = std(SFS2);    
            
            ROI_data = spm_detrend(ROI_data,1); %Detrend data
            
            %Filter ROI signal
            ROI_dataf = spm_filter(K,ROI_data);
            ROI_dataf = ROI_dataf-Xfilt*(pinv(Xfilt)*ROI_dataf);
          
            ROIvalidvoxels{roi} = size(ROI_dataf,2);
            ROIallvoxels{roi} = length(find(roiY==ROI_id));
            ROIrepre{roi} = mean(ROI_dataf,2);
            ROIcover{roi} = ROIvalidvoxels{roi} / ROIallvoxels{roi};
            
            %Calculate quality of representative signal as a mean correlation with all voxel in ROI
            ROI_CCx = corrcoef([ROIrepre{roi},ROI_dataf]);
            ROI_CC{roi} = mean(ROI_CCx(2:end,1));
            ROIs_info.repreCC{roi} = ROI_CC{roi};
            [repre_components reprevar] = PCA_MG(ROI_dataf);
            ROIs_info.repreV{roi} = reprevar(1);

            if ROIvalidvoxels{roi} > 1
                ROIs_info.repreVR{roi} = reprevar(1)/reprevar(2);
            else
                ROIs_info.repreVR{roi} = NaN;
            end
        %No valid voxels in ROI
        else 
           ROIrepre{roi} = NaN(nscans,1); 
           ROIvalidvoxels{roi} = 0;
           ROIallvoxels{roi} = length(find(roiY==ROI_id));
           ROIcover{roi} = 0;
           ROIs_info.SNR{roi} = NaN;
           ROIs_info.SNR2{roi} = NaN;
           ROIs_info.SNR3{roi} = NaN;
           ROIs_info.tSNRm{roi} = NaN;
           ROIs_info.tSNRs{roi} = NaN;
           ROIs_info.repreCC{roi} = NaN;
           ROIs_info.repreV{roi} = NaN;
           ROIs_info.repreVR{roi} = NaN;
        ROIs_info.SFS1m{roi} = NaN;
        ROIs_info.SFS1s{roi} = NaN;
        ROIs_info.SFS2m{roi} = NaN;
        ROIs_info.SFS2s{roi} = NaN;        
        end
        
        ROIs_info.ROIindex{roi} = roi;
        clear tSNR ROImean ROIstd ROI_data ROI_CCx ROI_data ROI_dataf
        
    end 

   RepreSigs = ROIrepre;
   ROIs_info.validvoxels = ROIvalidvoxels;
   ROIs_info.allvoxels = ROIallvoxels;
   ROIs_info.cover = ROIcover;
   ROIs_info.tSNR_GMm = tSNR_GMm;
   ROIs_info.tSNR_WMm = tSNR_WMm;
   ROIs_info.tSNR_GMs = tSNR_GMs;
   ROIs_info.tSNR_WMs = tSNR_WMs;
   ROIs_info.SNR_GM = SNR_GM;
   ROIs_info.SNR_WM = SNR_WM;
   ROIs_info.SNR_GM2 = SNR_GM2;
   ROIs_info.SNR_WM2 = SNR_WM2;
   ROIs_info.SNR_GM3 = SNR_GM3;
   ROIs_info.SNR_WM3 = SNR_WM3;
   ROIs_info.SFS_WM1m = SFS_WM1m;
   ROIs_info.SFS_WM1s = SFS_WM1s;
   ROIs_info.SFS_WM2m = SFS_WM2m;
   ROIs_info.SFS_WM2s = SFS_WM2s;
   ROIs_info.SFS_GM1m = SFS_GM1m;
   ROIs_info.SFS_GM1s = SFS_GM1s;
   ROIs_info.SFS_GM2m = SFS_GM2m;
   ROIs_info.SFS_GM2s = SFS_GM2s;
   ROIs_info.WMmean = mean(meanWM);
   ROIs_info.GMmean = mean(meanGM);
   ROIs_info.CSFmean = mean(meanCSF);
   ROIs_info.WMstd = mean(stdWM);
   ROIs_info.GMstd = mean(stdGM);
   ROIs_info.CSFstd = mean(stdCSF);
   ROIs_info.WMstdTS = mean(stdWMTS);
   ROIs_info.GMstdTS = mean(stdGMTS);
   ROIs_info.CSFstdTS = mean(stdCSFTS); 
   ROIs_info.CSFstdTSF = stdCSFTSF;
   ROIs_info.noisestdTSF = stdnoiseTSF;
   ROIs_info.noisemean = meannoise;
   ROIs_info.noisestd = stdnoise;
   ROIs_info.noisestdTS = stdnoiseTS;
   
   %Calculate SNS with respect to each individual noise voxel and save image with spatial distribution of this SNS calculation
   nrois = length(ROIs_info.cover);
   cover = cell2mat(ROIs_info.cover);
   vri = (cover>0.05); %Only valid indices of Repre Signals 5% of coverage is required here
   ValidRepreSigs = cell2mat(RepreSigs(vri)); %Prepare only valid signals
   vnrois = size(ValidRepreSigs,2);

   %Load definition of funstional networks within the ROI mask (prepared only for AAL)
   FCNdeffile = fullfile(pathstr,DefStruc.FuncNetworks_filename);
   FCNdef = load(FCNdeffile); 
   FCNdef = FCNdef.results;
      
   tmpROIcors = corrcoef(ValidRepreSigs);
   tempccind = find((triu(ones(vnrois))-eye(vnrois))>0); %Only upper part of matrix without diagonal will be stored
   ROIcors = tmpROIcors(tempccind);
      
   %Calculate ROIcors only within functional networs
   for net = 1:20
        netrois = find(FCNdef(:,net)>=0.5);
        fullvect = zeros(nrois,1);
        fullvect(netrois') = 1;
        validvect = fullvect(vri);
        netmatrix = validvect*validvect';
        netindices{net} = netmatrix(tempccind);
   end
   
   FCNcors=[];
   for net = 1:20
       FCNcors = [FCNcors; nonzeros(squeeze(abs(ROIcors.*netindices{net})))];
   end 
    
   tmpind = randperm(size(noise_dataf,2),min(1000,size(noise_dataf,2)));
   
   noise_dataf_reduced = noise_dataf(:,tmpind);
   SNS = zeros(1,size(noise_dataf_reduced,2));
   SNSt = zeros(1,size(noise_dataf_reduced,2));
   for ni = 1:size(noise_dataf_reduced,2)
      tmpmat =  [ValidRepreSigs noise_dataf_reduced(:,ni)];
      tmpcor = corrcoef(tmpmat);
      noisecors = tmpcor(end,1:end-1);
      nci(ni,:) = noisecors;
      SNSi = mean(abs(ROIcors))/mean(abs(noisecors));
      SNSi_FCN = mean(abs(FCNcors))/mean(abs(noisecors));
      [t,m1,m2,s1,s2,s,tp]=tt2calc(abs(ROIcors),abs(noisecors));
      SNSt(ni) = t;
      [t,m1,m2,s1,s2,s,tp]=tt2calc(abs(FCNcors),abs(noisecors));
      SNSt_FCN(ni) = t;
      SNS(ni) = SNSi;
      SNS_FCN(ni) = SNSi_FCN;
      
   end 

   SNSmean = mean(SNS(~isnan(SNS(:))));
   SNSstd = std(SNS(~isnan(SNS(:))));
   SNStmean = mean(SNSt(~isnan(SNSt(:))));
   SNStstd = std(SNSt(~isnan(SNSt(:))));
   SNSall= mean(abs(ROIcors))/mean(abs(nci(~isnan(nci(:)))));
   [t,m1,m2,s1,s2,s,tp]=tt2calc(abs(ROIcors),abs(nci(~isnan(nci(:)))));
   SNStall = t;
   
   SNSmean_FCN = mean(SNS_FCN(~isnan(SNS_FCN(:))));
   SNSstd_FCN = std(SNS_FCN(~isnan(SNS_FCN(:))));
   SNStmean_FCN = mean(SNSt_FCN(~isnan(SNSt_FCN(:))));
   SNStstd_FCN = std(SNSt_FCN(~isnan(SNSt_FCN(:))));
   SNSall_FCN= mean(abs(FCNcors))/mean(abs(nci(~isnan(nci(:)))));
   [t,m1,m2,s1,s2,s,tp]=tt2calc(abs(FCNcors),abs(nci(~isnan(nci(:)))));
   SNStall_FCN = t;
   
   ROIs_info.ROIcors = ROIcors;
   ROIs_info.ROIcors_mean = mean(ROIcors);
   ROIs_info.FCNcors = FCNcors;
   ROIs_info.FCNcors_mean = mean(FCNcors);
   ROIs_info.noisecors = noisecors;
   ROIs_info.noisecors_mean = mean(noisecors);
   ROIs_info.nci = nci;
   ROIs_info.SNS = SNS;
   ROIs_info.SNSmean = SNSmean;
   ROIs_info.SNSstd = SNSstd;
   ROIs_info.SNSall = SNSall;
   ROIs_info.SNSt = SNSt;
   ROIs_info.SNStmean = SNStmean;
   ROIs_info.SNStstd = SNStstd;
   ROIs_info.SNStall = SNStall;   
   
   ROIs_info.SNS_FCN = SNS_FCN;
   ROIs_info.SNSmean_FCN = SNSmean_FCN;
   ROIs_info.SNSstd_FCN = SNSstd_FCN;
   ROIs_info.SNSall_FCN = SNSall_FCN;
   ROIs_info.SNSt_FCN = SNSt_FCN;
   ROIs_info.SNStmean_FCN = SNStmean_FCN;
   ROIs_info.SNStstd_FCN = SNStstd_FCN;
   ROIs_info.SNStall_FCN = SNStall_FCN;    
end


% CoverThr specify the minimal cover of each ROI to be used within correlation matrix (to avoid regions withouth any signal or with very low SNR) - value between 0 and 1
function [CorrMat CorrMatFull] = PCorrelationMatrix(CoverThr, RepreSigs, ROIs_info, DefStruc)
    nrois = length(ROIs_info.cover);
    cover = cell2mat(ROIs_info.cover);
    vri = (cover>CoverThr); %Only valid indices of Repre Signals
    ValidRepreSigs = cell2mat(RepreSigs(vri)); %Prepare only valid signals
    vnrois = size(ValidRepreSigs,2);
    tempccind = find((triu(ones(vnrois))-eye(vnrois))>0); %Only upper part of matrix without diagonal will be stored
    
    cc = corrcoef(ValidRepreSigs);
    CorrMat.R = cc(tempccind);      
    CorrMat.Z = 0.5*log((1+CorrMat.R)./(1-CorrMat.R)); %Prepare Fisher z-transformed matrix too
    
    %Calculate full size correlation amtric (all ROIs) too
    AllRepreSigs = cell2mat(RepreSigs([1:nrois])); %Prepare all signals
    tempccind = find((triu(ones(nrois))-eye(nrois))>0);
    cc = corrcoef(AllRepreSigs);
    CorrMatFull.R = cc(tempccind);      
    CorrMatFull.Z = 0.5*log((1+CorrMatFull.R)./(1-CorrMatFull.R)); %Prepare Fisher z-transformed matrix too
end
    

function [outputdata] = mm_detrend(data,dim)    
    if dim~=1 
        data = data';
    end
    xdatamean = mean(data,1,'omitnan');
    xdatad = spm_detrend(data,1);
    outputdata = xdatad + xdatamean;
end   


function [components,var]=PCA_MG(Y)
% Function performs PCA analysis using svd. Based on SPM toolbox.
% Input:  Y - n signals for PCA decomposition - size Y:(signal length, n)
% Output: components - sorted N components - size components:(signal length, n)
%         var - explained variability [%] for every component - size components:(1, n)
% ------------------------------------------------------------------------
% Based on svd used in SPM toolbox
% Version 1.01 - use spm_detrend instead of ica_tb_detrend
% ------------------------------------------------------------------------
% Written by Martin Gajdos
% Last update 27.1.2016 14:36
% ------------------------------------------------------------------------

    [m,n]   = size(Y);

    if m > n
        [v s v] = svd(Y'*Y);
        s       = diag(s); %Eigenvalues
        
        s_perc=100*(s/sum(s));
        
        v       = v(:,1:numel(s_perc));
        u = NaN(m,numel(s_perc));
        for i=1:numel(s_perc)
            u(:,i)       = Y*v(:,i)/sqrt(s(i));
        end
        
    else
        [u s u] = svd(Y*Y');
        s       = diag(s); %Eigenvalues
        
        s_perc=100*(s/sum(s));
        
        u       = u(:,1:numel(s_perc));
        v = NaN(n,numel(s_perc));
        for i=1:numel(s_perc)
            v(:,i)       = Y'*u(:,i)/sqrt(s(i));
        end
    
    end

    d       = sign(sum(v));
    u       = u.*repmat(d,m,1);

    components=NaN(m,numel(s_perc));
    var=s_perc; %Explained variability [%]
    for i=1:numel(s_perc)
        components(:,i) = u(:,i)*sqrt(s(i)/n); %First eigenvariate {scaled - c.f. mean response}
    end
end

function [t,m1,m2,s1,s2,s,tp]=tt2calc(x,y)
    m1 = mean(x);
    m2 = mean(y);
    s1 = std(x);
    s2 = std(y);
    n1 = length(x);
    n2 = length(y);
    s = sqrt((((n1-1)*s1*s1)+((n2-1)*s2*s2))/(n1+n2-2));
    t = (m1-m2)/sqrt((s1*s1/n1)+(s2*s2/n2));
    tp = (m1-m2)/sqrt((s*s/n1)+(s*s/n2));
end 