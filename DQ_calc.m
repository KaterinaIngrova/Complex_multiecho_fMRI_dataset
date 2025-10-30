%% CALCULATION OF DATA QUALITY MAPS FOR fMRI DATA

function DQ_calc(fmri_data, DefStruc, normflag)

    % input:    fmri_data - NIFTI 4D fil with full path         
    %           DefStruc - other definitions
    %           normflag - 0 (default) if data are already normalized in MNI
    %                       space, 1 if data need to be normalized (it is
    %                       necessary for matching ROI/atlas template)

    if ~exist('normflag','var')
        normflag = 0;
    end
    
    
    if normflag == 1
        %normalize input data into MNI space
    end
    
    TR = DefStruc.TR;
    
    %Resample all masks into the same space according to first data file
    wmmask = DefStruc.WMmask_filename;
    gmmask = DefStruc.GMmask_filename;
    csfmask = DefStruc.CSFmask_filename;
    roimask = DefStruc.FuncNetMask_filename;
    datafiles = fmri_data;
    
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
    datav=spm_vol(datafiles);
    [dataY, ~]=spm_read_vols(datav);
    roiv=spm_vol(roimask);
    [roiY, ~]=spm_read_vols(roiv);
    
    
    %Find valid (inbrain) voxels and noise (background) voxel  
    nscans = size(datav,1);
    reduced_nscans = nscans; 
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
    valid_indices = (ibv_ii==min(nscans,reduced_nscans)); %Mask/indices of valid voxels 
    validindvec = find(valid_indices);
    noise_indices = (obv_ii==min(nscans,reduced_nscans)); %Mask/indices of noise voxels 
    noiseindvec = find(noise_indices);
    
    re_dataY = reshape(dataY,[],length(datav))'; %Reshape 4D data to 2D
    re_dataYd = mm_detrend(re_dataY,1);
    
    [DVARS, Stat] = DVARSCalc(re_dataY(:,validindvec)'); %Caluclate DVARS on inbrain voxels only
    ncrp_scans = find(Stat.pvals<0.05./(nscans-1) & Stat.DeltapDvar>5);
    DVARSout.DVARS = DVARS;
    DVARSout.Stat = Stat;
    DVARSout.ncrp_scans = ncrp_scans;
    
    re_roiY = reshape(roiY,[],1);
    
    WM_indices = (wmY>=DefStruc.WMmask_threshold)&valid_indices;
    CSF_indices = (csfY>=DefStruc.CSFmask_threshold)&valid_indices;
    GM_indices = (gmY>=DefStruc.GMmask_threshold)&valid_indices;
    WMindvec=find(WM_indices);
    GMindvec=find(GM_indices);
    CSFindvec=find(CSF_indices);
    
    re_tSNRmap = zeros(size(re_dataY,2),1);
    re_SFSmap = zeros(size(re_dataY,2),1);
    re_SFS2map = zeros(size(re_dataY,2),1);
    re_tSNRmapLFF = zeros(size(re_dataY,2),1);
    re_SFSmapLFF = zeros(size(re_dataY,2),1);    
    re_SNSmap = zeros(size(re_dataY,2),1);
    re_SNS2map = zeros(size(re_dataY,2),1);
    re_SNS3map = zeros(size(re_dataY,2),1);
    re_SNS3bmap = zeros(size(re_dataY,2),1);
    
    
    %Remove constant signal to avoid NaN in subsequent calculations
    temp_noise_data = re_dataY(:,noiseindvec);   
    ndstd = std(temp_noise_data);
    new_noiseindvec = noiseindvec(ndstd~=0);
    noiseindvec = new_noiseindvec;
    noise_data = temp_noise_data(:,ndstd~=0); %Final noise data
    
    noise_datad = mm_detrend(noise_data,1); %Remove linear trend
    noise_datadd = spm_detrend(noise_data,1); %Remove mean and linear trend
    meannoise = mean(mean(noise_datad,1));
    stdnoise = mean(std(noise_datad,0,2));
    stdnoiseTS = mean(std(noise_datad,0,1));    
    
    GloSIg = spm_global(datav);
    GloMean = mean(mean(re_dataY(:,validindvec)));
    
    Faxis = (0:nscans-1)/nscans/TR;
    complete_range = Faxis > 0 & Faxis < (1/TR/2);
    bold_range = Faxis > 0.01 & Faxis < 0.2 & Faxis < (1/TR/2);
    noise_range = logical(complete_range - bold_range);
    PSdata = 2*abs(fft(re_dataYd))/nscans; %Power spectrum of data
    PSdata(1,:) = 0; %Discard DC
    completeP = sum(PSdata(complete_range,:));
    boldP = sum(PSdata(bold_range,:));
    noiseP = sum(PSdata(noise_range,:));
    
    re_tSNRmap(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan') ./ std(re_dataYd(:,validindvec),0,1,'omitnan')) * sqrt(length(datav)) ;
    re_tSNRmapLFF(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan')./GloMean) ./ noiseP(:,validindvec) * sqrt(length(datav)) ;
    
    CSF_data = re_dataY(:,CSFindvec);
    CSF_datad = mm_detrend(CSF_data,1); %Remove linear drift
    meanCSF = mean(mean(CSF_datad,1));
    stdCSF = mean(std(CSF_datad,0,2));
    stdCSFTS = mean(std(CSF_datad,0,1));    
    WM_data = re_dataY(:,WMindvec);
    WM_datad = mm_detrend(WM_data,1); %Remove linear drift
    meanWM = mean(mean(WM_datad,1));
    stdWM = mean(std(WM_datad,0,2));
    stdWMTS = mean(std(WM_datad,0,1));
    GM_data = re_dataY(:,GMindvec);
    GM_datad = mm_detrend(GM_data,1);
    meanGM = mean(mean(GM_datad,1));
    stdGM = mean(std(GM_datad,0,2));
    stdGMTS = mean(std(GM_datad,0,1));
    
    [WMcomp,WMvar]=PCA_MG(WM_datad);
    [CSFcomp,CSFvar]=PCA_MG(CSF_datad);
    
    X=[WMcomp(:,1:5) CSFcomp(:,1:5) ones(nscans,1)];
    cw = [1 1 1 1 1 1 1 1 1 1 0];
    betas=(pinv(X)*re_dataYd);
    re_dataYdf = re_dataYd-X(:,1:10)*betas(1:10,:);
    
    betasn=(pinv(X)*noise_datadd);
    noise_dataddf = noise_datadd-X(:,1:10)*betasn(1:10,:);
    
    PSdataf = 2*abs(fft(re_dataYdf))/nscans; %Power spectrum of data
    PSdataf(1,:) = 0; %Discard DC
    completePf = sum(PSdataf(complete_range,:));
    boldPf = sum(PSdataf(bold_range,:));
    noisePf = sum(PSdataf(noise_range,:));
    
    re_tSNRmapLFF(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan')./GloMean) ./ (noiseP(:,validindvec)+(boldP(:,validindvec)-boldPf(:,validindvec)));
    
    
    re_SFSmap(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan')./GloMean).*(std(re_dataYd(:,validindvec),0,1,'omitnan')./stdnoiseTS);
    re_SFSmapLFF(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan')./GloMean).*boldP(:,validindvec)./(noiseP(:,validindvec)+(boldP(:,validindvec)-boldPf(:,validindvec)));
    re_SFS2map(validindvec) = (mean(re_dataYd(:,validindvec),1,'omitnan')./GloMean).*(std(re_dataYd(:,validindvec),0,1,'omitnan')./stdCSFTS);
    
    
    %Calculate roi representatives for within-brain correclations
    roivalues = unique(re_roiY(re_roiY~=0 & ~isnan(re_roiY)));
    nrois = size(roivalues,1);
    
    roi_repre = [];
    roi_repref = [];
    roi_voxcount=[];
    roi_mean=[];
    roi_meanf=[];
    for iroi = 1:nrois
        roi_indices = re_roiY==roivalues(iroi);
        roi_data = re_dataYd(:,roi_indices);
        roi_dataf = re_dataYdf(:,roi_indices);
        if ~isempty(roi_data)
            roi_voxcount(iroi) = size(roi_data,2);
            roi_datad = spm_detrend(roi_data,1);
            roi_datadf = spm_detrend(roi_dataf,1);
            roi_mean(iroi) = mean(mean(roi_data,1,'omitnan'));
            tmprepre = mean(roi_datad,2,'omitnan');
            tmprepref = mean(roi_datadf,2,'omitnan');
            roi_repre = [roi_repre tmprepre];
            roi_repref = [roi_repref tmprepref];
        end
    end
    
    %Select noise_repre from all noise voxels (for faster computation of correlations) 
    noiseid = floor(rand(1,min(1000,size(noise_datadd,2)))*(size(noise_datadd,2)-1))+1;
    noise_repre = noise_datadd(:,noiseid);
    noise_repref = noise_dataddf(:,noiseid);
    
    %Faster version with matrix calculation of pearson correlation
    X = re_dataYd(:,validindvec); %Only valid data
    Y = roi_repre;
    
    X = bsxfun(@minus,X,nansum(X,1)./size(X,1)); 
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    
    %Normalize by the L2-norm (Euclidean) of Rows:
    X = X.*repmat(sqrt(1./max(eps,nansum(abs(X).^2,1))),[size(X,1),1]); 
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    %Compute Pair-wise Correlation Coefficients:
    r = (X'*Y);
    reprecor = mean(abs(r),2);
    repredev = std(abs(r),0,2);
    
    Y = noise_repre;
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    r = (X'*Y);
    noisecor = mean(abs(r),2);
    noisedev = std(abs(r),0,2);
    
    Y = noise_repref;
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    r = (X'*Y);
    noisecorf = mean(abs(r),2);
    noisedevf = std(abs(r),0,2);
    
    sns2 = reprecor./noisecor;
    sns = sns2 .* (noisedev./repredev);
    re_SNSmap(validindvec) = sns2;
    re_SNS2map(validindvec) = sns;
    
    %Repeat for filtered data - just inbrain correalations
    X = re_dataYdf(:,validindvec); %Only valid data
    Y = roi_repref;
    
    X = bsxfun(@minus,X,nansum(X,1)./size(X,1)); 
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    
    %Normalize by the L2-norm (Euclidean) of Rows
    X = X.*repmat(sqrt(1./max(eps,nansum(abs(X).^2,1))),[size(X,1),1]); 
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    %Compute Pair-wise Correlation Coefficients
    r = (X'*Y);
    reprecorf = mean(abs(r),2);
    repredevf = std(abs(r),0,2);   
    
    deltacor = (reprecor-reprecorf)./reprecor;
    SNScorfact = 1-deltacor;
    re_SNS3map(validindvec) = (reprecor.*SNScorfact)./noisecor;
    re_SNS3bmap(validindvec) = reprecorf./noisecorf;
    
    SNSmap = reshape(re_SNSmap,datav(1).dim); %reshape 1D data back to 3D
    SNS2map = reshape(re_SNS2map,datav(1).dim);
    SNS3map = reshape(re_SNS3map,datav(1).dim);
    SNS3bmap = reshape(re_SNS3bmap,datav(1).dim);
    tSNRmap = reshape(re_tSNRmap,datav(1).dim);
    SFSmap = reshape(re_SFSmap,datav(1).dim);
    tSNRmapLFF = reshape(re_tSNRmapLFF,datav(1).dim);
    SFSmapLFF = reshape(re_SFSmapLFF,datav(1).dim);
    SFS2map = reshape(re_SFS2map,datav(1).dim);
    
    %save DQ maps
    [olddatadir, datafn , ~] = fileparts(datav(1).fname);
    datadir=fullfile(olddatadir,'DQoutputs');
    if ~exist(datadir,'dir')
        mkdir(datadir);
    end
    cd(datadir)
    [fdir,fname,fext]=fileparts(fmri_data);
    save(['DVARSout' fname '.mat'],'DVARSout');
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['tSNRmap_' datafn '.nii']);
    outputv.descrip=['tSNR map (tSNR multiplied with sqrt(nscans)'];
    outputv = spm_write_vol(outputv,tSNRmap);   
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['tSNRLFFmap_' datafn '.nii']);
    outputv.descrip=['tSNR_LFF map (multiplied with sqrt(nscans)'];
    outputv = spm_write_vol(outputv,tSNRmapLFF);      
    
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SFSmap_' datafn '.nii']);
    outputv.descrip=['SFS map based on background noise'];
    outputv = spm_write_vol(outputv,SFSmap);  
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SFSLFFmap_' datafn '.nii']);
    outputv.descrip=['SFS map based on LFF bold2noise ratio'];
    outputv = spm_write_vol(outputv,SFSmapLFF);     
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SFS2map_' datafn '.nii']);
    outputv.descrip=['SFS map based on CSF noise'];
    outputv = spm_write_vol(outputv,SFS2map);  
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SNSmap_' datafn '.nii']);
    outputv.descrip=['SNS map based on ratio'];
    outputv = spm_write_vol(outputv,SNSmap);  
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SNS2map_' datafn '.nii']);
    outputv.descrip=['SNS map based on ratio normalized by standard deviation'];
    outputv = spm_write_vol(outputv,SNS2map);  
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SNS3map_' datafn '.nii']);
    outputv.descrip=['SNS map based on ratio corrected by nuisance effects'];
    outputv = spm_write_vol(outputv,SNS3map); 
    
    outputv = datav(1);
    outputv.dt = [spm_type('float32') spm_platform('bigend')];
    outputv.fname = fullfile(datadir,['SNS3bmap_' datafn '.nii']);
    outputv.descrip=['SNS map based on ratio corrected by nuisance effects'];
    outputv = spm_write_vol(outputv,SNS3bmap); 
    
    end
    
    function [r,p] = fast_corr(X,Y)
    % [R,P] = FAST_CORR(X,Y)
    % 
    % Enables the quick vectorized computation of pair-wise correlations 
    % between corresponding columns of two large matrices.
    % 
    % -Around 6.75 times faster than using corr.m in a for loop for 100 
    %   observations, 1000 variables
    % -Over 11 times faster than computing the full correlation matrix using
    %   corr.m
    % 
    % Inputs:
    % X & Y,    size [n_observations x n_variables] matrices/vectors; both must
    %           be of equal size to enable computation of pair-wise
    %           correlations, column-by-column
    % 
    % Outputs:
    % r,        a vector of Pearson product-moment correlation coefficients of
    %           length equal to the number of columns in X and Y; columns of X 
    %           and Y are correlated pair-wise
    % p,        a vector of p-values corresponding to r
    % 
    % Author:  Elliot Layden, The University of Chicago, 2017
    
    %Check Data Sizes:    [r1,c1] = size(X); [r2,c2] = size(Y); 
    if r1~=r2 
        error('''X'' and ''Y'' must contain the same number of rows/observations.')
    end
    if c1~=c2 
        error('''X'' and ''Y'' must contain the same number of columns/variables.')
    end
    
    %Pair-wise Removal of Rows w/ NaN's:
    %   Note: this removes any NaN-containing rows from entire matrices;
    %       otherwise speed would be sacrificed
    if any(isnan(X(:))) || any(isnan(Y(:)))
        nan1 = isnan(X); nan2 = isnan(Y);
        nans = sum([nan1,nan2],2)>0;
        X(nans) = []; Y(nans) = [];
    end
    
    %De-mean Columns:
    X = bsxfun(@minus,X,nansum(X,1)./size(X,1)); 
    Y = bsxfun(@minus,Y,nansum(Y,1)./size(Y,1));
    
    %Normalize by the L2-norm (Euclidean) of Rows:
    X = X.*repmat(sqrt(1./max(eps,nansum(abs(X).^2,1))),[size(X,1),1]); 
    Y = Y.*repmat(sqrt(1./max(eps,nansum(abs(Y).^2,1))),[size(Y,1),1]);
    
    %Compute Pair-wise Correlation Coefficients:
    r = nansum(X.*Y);
    
    %Calculate p-values if requested:
    if nargout==2
        t = (r.*sqrt(r1-2))./sqrt(1-r.^2);
        p = 2*tcdf(abs(t),(r1-2),'upper');
    end
    
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
        s = diag(s); %Eigenvalues
        
        s_perc=100*(s/sum(s));
        
        v = v(:,1:numel(s_perc));
        u = NaN(m,numel(s_perc));
        for i=1:numel(s_perc)
            u(:,i) = Y*v(:,i)/sqrt(s(i));
        end
    
    else
        [u s u] = svd(Y*Y');
        s       = diag(s); % s - eigenvalues
        
        s_perc=100*(s/sum(s));
        
        u       = u(:,1:numel(s_perc));
        v = NaN(n,numel(s_perc));
        for i=1:numel(s_perc)
            v(:,i) = Y'*u(:,i)/sqrt(s(i));
        end
    end
    
    d = sign(sum(v));
    u = u.*repmat(d,m,1);
    
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