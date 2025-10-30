%% PROCESS MULTI ECHO DATA
function Process_MultiEcho

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

    %By number 1 specify what preprocesing and calculation should be done
    ToDo.realign = 1;
    ToDo.unwarp = 0;  
    ToDo.calcsumimage = 0;
    ToDo.calcSNRcomp = 0;
    ToDo.calctSNRcomp = 1;

    %Thresholds foe finding outbrain ans inbrain voxels
    outbrain_thr = 0.20; %Background
    inbrain_thr = 0.25; %Valid data for analysis


    spm_defaults;
    global defaults;

    %For loop across subjects
    for sub=1:nsub 

        %For loop across sessions
        for ses=1:nses
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
            % Realign
            %-----------------------------------------------------------------
            if ToDo.realign == 1      
                P     = [];
                selected_scans = [];
                direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func']; %Path to functional data   
                q   = spm_select('FPList',direc{ses},['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-2_part-mag_bold\.nii$']); %Selecting functional data of second echo

                P   = strvcat(P,q);
                for i=1:size(P,1)
                    selected_scans{i,1} = [P(i,:)];
                end

               %Decide for realign step with unwarping
               if ToDo.unwarp == 1 %Do realing and unwarp
                    ReaUnwPrep = 'u';
                    for ee=1:NoOfEchoes   
                        RUW.data(ee).scans=selected_scans(ee);
                        RUW.data(ee).pmscan='';
                    end

                %Set parameters of realign
                RUW.eoptions   = struct(  'quality',     0.9,...
                                     'sep',         4,...
                                     'fwhm',        5,...
                                     'rtm',         0,...
                                     'einterp',     2,...
                                     'ewrap',       [0 0 0],...
                                     'weight',      '');

                RUW.uweoptions   = struct(  'basfcn',     [12 12],...
                                     'regorder',         1,...
                                     'lambda',        100000,...
                                     'jm',          0,...
                                     'fot',         [4 5],...
                                     'sot',         [],...
                                     'uwfwhm',      4,...
                                     'rem',         1,...
                                     'noi',         5,...
                                     'expround',    'Average');


                RUW.uwroptions = struct(  'uwwhich',     [2 1],...
                                     'rinterp',     4,...
                                     'wrap',        [0 0 0],...
                                     'mask',        1,...
                                     'prefix',      'u');

                matlabbatch{1}.spm.spatial.realignunwarp = RUW;                      
                spm_jobman('run',matlabbatch);    
                clear matlabbatch RUW P q selected_scans dfiles

            %Do only realign without unwarping 
            else 
                ReaUnwPrep = 'r';
                for ee=1:1   
                    dfiles(ee)={selected_scans(ee)};
                end
                REW.data = dfiles;

                REW.eoptions   = struct(  'quality',     0.9,...
                                     'sep',         4,...
                                     'fwhm',        5,...
                                     'rtm',         0,...
                                     'interp',     2,...
                                     'wrap',       [0 0 0],...
                                     'weight',      '');


                REW.roptions = struct(  'which',     [2 1],...
                                     'interp',     4,...
                                     'wrap',        [0 0 0],...
                                     'mask',        1,...
                                     'prefix',      'r');

                matlabbatch{1}.spm.spatial.realign.estwrite = REW;                      
                spm_jobman('run',matlabbatch);    
                clear matlabbatch REW P q selected_scans dfiles  


                %Realign of other echoes and phase images according to realign of the second echo        
                direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func'];
                q   = spm_select('FPList',direc{ses},['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-2_part-mag_bold.nii$']);


                for e = [1 3] 
                    otherdirec{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func'];  
                    otherq   = spm_select('FPList',otherdirec{ses},['^' subjects{sub,2} '.*' sessions{ses,1} '_.*echo-' num2str(e) '_part-mag_bold.nii$']);
                    rewrite_mat_reslice(q,otherq)     
                    clear otherdirec otherq

                end     

                if PHflag 
                    for e = [1 2 3] 
                        otherdirec{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func'];
                        otherq   = spm_select('FPList',otherdirec{ses},['^' subjects{sub,2}  '.*' sessions{ses,1} '_.*echo-' num2str(e) '_part-phase_bold.nii$']);
                        rewrite_mat_reslice(q,otherq)     
                        clear otherdirec otherq
                    end  
                end
            end
        end



            %-----------------------------------------------------------------     
            % Combining echoes
            %-----------------------------------------------------------------
            %Load realigned files
            direc{ses}  = [cwd pathdelim subjects{sub,1} pathdelim 'func'];     
            q   = spm_select('FPList',direc{ses},['^r.*' sessions{ses,1} '_echo-._part-mag_bold.nii$']);        
            fv = spm_vol(q);
            nfiles = size(fv,1);

            nscans = nfiles/NoOfEchoes;

            dim = fv(1).dim;

            meY=cell(3);
            for echo=1:NoOfEchoes 
                fileindices = ((echo-1)*nscans)+1:(echo*nscans);
                [meY{echo}, ~] = spm_read_vols(fv(fileindices));
            end

            %Find outbrain and inbrain voxels in all data (all scans and echoes)
            obv_ii=zeros(size(meY{1}(:,:,:,1)));
            ibv_ii1e=zeros(size(meY{1}(:,:,:,1)));
            ibv_ii=zeros(size(meY{1}(:,:,:,1)));
            for scan=1:min(nscans,100)
                for echo=1:NoOfEchoes 
                    Yi = meY{echo}(:,:,:,scan);
                    meanY = mean(Yi(~isnan(Yi)));
                    othr=meanY*outbrain_thr;
                    ithr=meanY*inbrain_thr;
                    obv_ii = obv_ii + ((Yi<othr)&(~isnan(Yi)));
                    ibv_ii = ibv_ii + ((Yi>ithr)&(~isnan(Yi)));
                    if echo==1
                        ibv_ii1e = ibv_ii1e + ((Yi>ithr)&(~isnan(Yi)));
                    end
                end
            end
            obv_indices = (obv_ii==round(min(nscans,100)*NoOfEchoes)); 
            ibvAe_indices = (ibv_ii==round(min(nscans,100)*NoOfEchoes)); 
            ibv1e_indices = (ibv_ii1e==round(min(nscans,100)));

            ibv_indices = ibv1e_indices;

            nvox = dim(1)*dim(2)*dim(3);

            %Reshape data for faster calculations
            re_ibv_indices = reshape(ibv_indices,nvox,1);
            re_obv_indices = reshape(obv_indices,nvox,1);
            for echo=1:NoOfEchoes 
                re_meY{echo} = reshape(meY{echo},[],nscans);
            end

            %Calculate global outbrain signal variations
            for scan=1:nscans
                for echo=1:NoOfEchoes 
                    Yi = meY{echo}(:,:,:,scan);
                    bckg_sig(scan,echo) = mean(Yi(obv_indices));
                end
            end
            %Calculate reference signal value from each echo
            for echo=1:NoOfEchoes 
                ref_sig(echo) = mean(bckg_sig(:,echo));
            end


           %Prepare variables for composite echo
           if ToDo.calcsumimage == 1
               sumY=zeros(dim(1),dim(2),dim(3),nscans);
               re_sumY=reshape(sumY,[],nscans);
           end

           if ToDo.calctSNRcomp == 1
               tsnrY=zeros(dim(1),dim(2),dim(3),nscans);
               re_tsnrY=reshape(tsnrY,[],nscans);
           end

           if ToDo.calcSNRcomp == 1
               snrY=zeros(dim(1),dim(2),dim(3),nscans);
               re_snrY=reshape(snrY,[],nscans);
           end

           nvox = dim(1)*dim(2)*dim(3);
           nsegments = dim(3); 
           segsize = ceil(nvox/nsegments);
           for seg=1:nsegments
                segindices{seg} = ((seg-1)*segsize)+1:min((seg*segsize),nvox); 
           end

           if ToDo.calctSNRcomp == 1 
           %Calculate tSNR for each voxel   
                sumtSNR=zeros(nvox,1);
                for echo=1:NoOfEchoes
                    for seg=1:nsegments
                        re_tSNR{echo}(segindices{seg},1) = mean(re_meY{echo}(segindices{seg},:),2,'omitnan')./std(re_meY{echo}(segindices{seg},:),0,2,'omitnan');
                        sumtSNR(segindices{seg},1) = sumtSNR(segindices{seg},1) + (re_tSNR{echo}(segindices{seg},1)*TEs(echo));
                    end
                end                    
           end

           %Create composite multiecho image using sumation or average                 
           if ToDo.calcsumimage == 1
               sumY=zeros(size(re_meY{1}));
                %First calculate simple sum or mean from all echoes
                for echo=1:NoOfEchoes
                    for seg=1:nsegments
                        sumY(segindices{seg})=sumY(segindices{seg})+re_meY{echo}(segindices{seg});
                    end
                end
           end

           if ToDo.calctSNRcomp == 1
                %Create weighted mean using tSNR
                re_tsnrY=zeros(size(re_meY{1}));
                for seg=1:nsegments
                    for echo=1:NoOfEchoes
                        re_tsnrY(segindices{seg},:) = re_tsnrY(segindices{seg},:)+(re_meY{echo}(segindices{seg},:).*repmat(re_tSNR{echo}(segindices{seg},1)*TEs(echo),[1 nscans]));
                    end
                    re_tsnrY(segindices{seg},:) = re_tsnrY(segindices{seg},:)./repmat(sumtSNR(segindices{seg},1),[1 nscans]);
                end
           end

            %Write imporatant variables    
            [pathstr,name,ext]=fileparts(fv(1).fname);
            varfname=fullfile(pathstr,'ProcessMultiEcho_vars.mat');
            save(varfname,'NoOfEchoes','nscans','nfiles','fv','q','dim');

               if ToDo.calctSNRcomp == 1
                   for echo=1:NoOfEchoes
                        tSNR{echo} = reshape(re_tSNR{echo},dim(1),dim(2),dim(3));
                   end
                   save(varfname,'tSNR','sumtSNR','-append');
               end

               if ToDo.calcSNRcomp == 1
                   for echo=1:NoOfEchoes
                        SNR{echo} = reshape(re_SNR{echo},dim(1),dim(2),dim(3));
                   end
                   save(varfname,'SNR','sumSNR','ref_sig','-append');
               end

            %Write output files 
            if ToDo.calctSNRcomp == 1
                %Write tSNR files for each echo
                for echo=1:NoOfEchoes   
                    Vo=fv(1);
                    fname=Vo.fname;
                    [pathstr,name,ext]=fileparts(fname);
                    Vo.fname=fullfile(pathstr,['tSNR_' num2str(echo) '_' name ext]);
                    Vo.descrip=['tSNR_' num2str(echo)];
                    Vo=spm_write_vol(Vo,tSNR{echo});
                end
            end

             if ToDo.calcSNRcomp == 1
                %Write SNR files for each echo
                for echo=1:NoOfEchoes   
                    Vo=fv(1);
                    fname=Vo.fname;
                    [pathstr,name,ext]=fileparts(fname);
                    Vo.fname=fullfile(pathstr,['SNR_' num2str(echo) '_' name ext]);
                    Vo.descrip=['SNR_' num2str(echo)];
                    Vo=spm_write_vol(Vo,SNR{echo});
                end    
             end


            %Write 4D composite ME data 
            if ToDo.calcsumimage == 1
                sumY = reshape(re_sumY,dim(1),dim(2),dim(3),nscans); %Sumation of all echoes
                data=sumY;
                for ii=1:nscans
                    V(ii)=fv(ii);
                    fname=V(ii).fname;
                    [pathstr,name,ext]=fileparts(fname);
                    V(ii).fname=fullfile(pathstr,['sue_' name ext]);
                    V(ii).dim=size(data(:,:,:,ii));
                    V(ii).descrip='composite multiecho - sum';
                    spm_write_vol(V(ii),data(:,:,:,ii));
                end
            end


            if ToDo.calcSNRcomp == 1 
                snrY = reshape(re_snrY,dim(1),dim(2),dim(3),nscans); %Weighted average according to SNR
                data=snrY;
                for ii=1:nscans
                    V(ii)=fv(ii);
                    fname=V(ii).fname;
                    [pathstr,name,ext]=fileparts(fname);
                    V(ii).fname=fullfile(pathstr,['sne_' name ext]);
                    V(ii).dim=size(data(:,:,:,ii));
                    V(ii).descrip='composite multiecho - SNR weighted';
                    spm_write_vol(V(ii),data(:,:,:,ii));
                end
            end

            if ToDo.calctSNRcomp == 1
                tsnrY = reshape(re_tsnrY,dim(1),dim(2),dim(3),nscans); %Weighted average according to tSNR
                data=tsnrY;
                for ii=1:nscans
                    V(ii)=fv(ii);
                    fname=V(ii).fname;
                    [pathstr,name,ext]=fileparts(fname);
                    V(ii).fname=fullfile(pathstr,['tse_' name ext]);
                    V(ii).dim=size(data(:,:,:,ii));
                    V(ii).descrip='composite multiecho - tSNR weighted';
                    spm_write_vol(V(ii),data(:,:,:,ii));
                end
            end

        clear  meY snrY tsnrY sumY tSNR SNR re_meY re_tsnrY;
       end 
    end
end 


function [dataout] = removeoutliers(datain)
%REMOVEOUTLIERS   Remove outliers from data using the Thompson Tau method.
%   For vectors, REMOVEOUTLIERS(datain) removes the elements in datain that
%   are considered outliers as defined by the Thompson Tau method. This
%   applies to any data vector greater than three elements in length, with
%   no upper limit (other than that of the machine running the script).
%   Additionally, the output vector is sorted in ascending order.
%
%   Example: If datain = [1 34 35 35 33 34 37 38 35 35 36 150]
%
%   then removeoutliers(datain) will return the vector:
%       dataout = 33 34 34 35 35 35 35 36 37 38
%
%   See also MEDIAN, STD, MIN, MAX, VAR, COV, MODE.
%   This function was written by Vince Petaccio on July 30, 2009.
n=length(datain); %Determine the number of samples in datain
if n < 3
    display(['ERROR: There must be at least 3 samples in the' ...
        ' data set in order to use the removeoutliers function.']);
else
    S=std(datain); %Calculate S, the sample standard deviation
    xbar=mean(datain); %Calculate the sample mean
    %tau is a vector containing values for Thompson's Tau
    tau = [1.150 1.393 1.572 1.656 1.711 1.749 1.777 1.798 1.815 1.829 ...
        1.840 1.849 1.858 1.865 1.871 1.876 1.881 1.885 1.889 1.893 ...
        1.896 1.899 1.902 1.904 1.906 1.908 1.910 1.911 1.913 1.914 ...
        1.916 1.917 1.919 1.920 1.921 1.922 1.923 1.924];
    %Determine the value of S times Tau
    if n > length(tau)
        TS=1.960*S; %For n > 40
    else
        TS=tau(n)*S; %For samples of size 3 < n < 40
    end
    %Sort the input data vector so that removing the extreme values
    %becomes an arbitrary task
    dataout = sort(datain);
    %Compare the values of extreme high data points to TS
    while abs((max(dataout)-xbar)) > TS 
        dataout=dataout(1:(length(dataout)-1));
        %Determine the NEW value of S times Tau
        S=std(dataout);
        xbar=mean(dataout);
        if length(dataout) > length(tau)
            TS=1.960*S; %For n > 40
        else
            TS=tau(length(dataout))*S; %For samples of size 3 < n < 40
        end
    end
    %Compare the values of extreme low data points to TS.
    %Begin by determining the NEW value of S times Tau
        S=std(dataout);
        xbar=mean(dataout);
        if length(dataout) > length(tau)
            TS=1.960*S; %For n > 40
        else
            TS=tau(length(dataout))*S; %For samples of size 3 < n < 40
        end
    while abs((min(dataout)-xbar)) > TS 
        dataout=dataout(2:(length(dataout)));
        %Determine the NEW value of S times Tau
        S=std(dataout);
        xbar=mean(dataout);
        if length(dataout) > length(tau)
            TS=1.960*S; %For n > 40
        else
            TS=tau(length(dataout))*S; %For samples of size 3 < n < 40
        end
    end
end
end



function [Y] = get_first_comp(y)
%Compute regional response in terms of first eigenvariate
%-----------------------------------------------------------------------
[m n]   = size(y);
if m > n
	[v s v] = svd(y'*y);
	s       = diag(s);
	v       = v(:,1);
	u       = y*v/sqrt(s(1));
else
	[u s u] = svd(y*y');
	s       = diag(s);
	u       = u(:,1);
	v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y      = u*sqrt(s(1)/n);
end

function rewrite_mat_reslice(files,otherfiles)
        mujfile_mask = spm_vol(files);

        mujfile_other = spm_vol(otherfiles);
        for scan = 1:size(mujfile_other,1)
        
            mujfile_other(scan).mat= mujfile_mask(scan).mat;
     
            
            %Update of transformation matrix in images (write the new matrix)
            spm_get_space([mujfile_other(scan).fname ',' num2str(mujfile_other(scan).n)], mujfile_mask(scan).mat);
        
            
        end
        
    matlabbatch{1}.spm.spatial.realign.write.data = {otherfiles};
    matlabbatch{1}.spm.spatial.realign.write.roptions.which = [2 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
    clear matlabbatch;      
        
end