% 2026 TS
% This script fully anonymizes BIDS dataset, developed and tailored to
% MAFIL data (CEITEC). Currently it supports MR, physio, stim.logs, other
% modalities not fully implemented.
% Prerequisities:
% 1) SPM12 toolbox in Matlab path, tested with v7771 - https://www.fil.ion.ucl.ac.uk/spm/
% LUT table to map OLD to NEW IDs has to be specified below, to be fully
% anonymized, the LUT table has to be deleted at the end of the
% anonymization process.

pth='e:\PROJECT\bids\';

todoLUT={          % use hash, randn to generace, than recode
    '1234A'	'001'
    };

% folder = uigetdir(pwd, 'Select a BIDS folder to ANONYMIZE');
% if folder==0, return; end
% cd(folder);

% Find files to be dealt with
todofiles=dir(fullfile(pth,'**/*.*'));
todofiles([todofiles.isdir]) = []; % skip folders now, deal with them later
fldsanon={
    'PatientName','PatientID','PatientWeight','PatientBirthDate','SeriesInstanceUID','StudyInstanceUID','AcquisitionDateTime','ReferringPhysicianName','InstitutionName','ProcedureStepDescription'
    }; % Field to be removed from MR metadata JSONs


% FIRST RUN, check file contents, rename files later in second run, rename folders in third run
h = waitbar(0,'Processing');
for ind=1:size(todofiles,1) 
    waitbar(ind/size(todofiles,1))
    
    % Warning, use manual check for those filetypes
    if any(regexpi(todofiles(ind).name,'^task-.*_bold.json$|^.bidsignore$|^dataset_description.json$|^participants.json$','once'))        % .bidsignore
        fprintf(1,'Please check file "%s" manually for possible sensitive/personal information\n',todofiles(ind).name)
        continue 
       
    elseif any(regexpi(todofiles(ind).name,'^.*.(bvec|bval)$'))                                 % skip those files
        continue
       
    elseif any(regexpi(todofiles(ind).name,'^participants.tsv$'))                               % participants.tsv
             parttsv=spm_load(fullfile(todofiles(ind).folder,todofiles(ind).name));
             for indd=1:numel(parttsv.participant_id)
                 for ind2=1:size(todoLUT,1)
                     if any(regexp(parttsv.participant_id{indd},todoLUT{ind2,1}))  % match participant name/row
                         parttsv.participant_id{indd}=regexprep(parttsv.participant_id{indd},todoLUT{ind2,1},todoLUT{ind2,2}); % replace string based on LUT
                         break
                     end
                 end
             end
             spm_save(fullfile(todofiles(ind).folder,todofiles(ind).name),parttsv)
             continue
             
    elseif any(regexpi(todofiles(ind).name,'^_list_of_original_sequences.tsv$'))                 % _list_of_original_sequences.tsv
             listorigseqtsv=spm_load(fullfile(todofiles(ind).folder,todofiles(ind).name));
             for indd=1:numel(listorigseqtsv.participant_id)
                 for ind2=1:size(todoLUT,1)
                     if any(regexp(listorigseqtsv.participant_id{indd},todoLUT{ind2,1}))  % match participant name/row
                         listorigseqtsv.participant_id{indd}=regexprep(listorigseqtsv.participant_id{indd},todoLUT{ind2,1},todoLUT{ind2,2}); % replace string based on LUT
                         listorigseqtsv.BIDS_rename{indd}=regexprep(listorigseqtsv.BIDS_rename{indd},todoLUT{ind2,1},todoLUT{ind2,2}); % replace string based on LUT
                         break
                     end
                 end
             end
             spm_save(fullfile(todofiles(ind).folder,todofiles(ind).name),listorigseqtsv)
             continue
             
    
    elseif any(regexpi(todofiles(ind).name,'(?!.*(physio|events).*)^sub.*.json$'))                      %  JSON (except PHYSIO JSON)
         for ind2=1:size(todoLUT,1)
             if any(regexp(todofiles(ind).name,todoLUT{ind2,1}))
                ofnfp=fullfile(todofiles(ind).folder,todofiles(ind).name); % old filename full path
                aa=spm_jsonread(ofnfp);
                for ine=1:size(fldsanon,2) % anonymization by removing tags from json
                    if isfield(aa,fldsanon{ine})
                        aa=rmfield(aa,fldsanon{ine});
                    end
                end
                aa.AccessionNumber=regexprep(aa.AccessionNumber,todoLUT{ind2,1},todoLUT{ind2,2});

                % IntendedFor
                if isfield(aa, 'IntendedFor')
                    aa.IntendedFor = regexprep(aa.IntendedFor, todoLUT{ind2,1}, todoLUT{ind2,2});
                end


                spm_jsonwrite(ofnfp,aa,struct('indent',' ')) % overwrite old, rename to new filename later
                clear ofpfn aa
                break
             end
         end
         continue
         
    elseif any(regexpi(todofiles(ind).name,'^.*.nii$'))                       %  NII
        for ind2=1:size(todoLUT,1)
            if any(regexp(todofiles(ind).name,todoLUT{ind2,1}))
                ofnfp=fullfile(todofiles(ind).folder,todofiles(ind).name); % old filename full path
                nii = nii_tool('load',ofnfp); % load the file
                if ~isempty(nii.hdr.db_name)
                    nii.hdr.db_name=regexprep(nii.hdr.db_name,todoLUT{ind2,1},todoLUT{ind2,2});
                    nii_tool('save', nii); % save it, overwrite old, rename to new filename later
                end
                clear ofpfn nii
                break
            end
        end
         continue
         
    elseif any(regexpi(todofiles(ind).name,'^.*.mat$'))                          %  MAT-files
        fprintf(1,'MAT-file "%s" detected, manually search for SENSITIVE data\n',todofiles(ind).name)
        ischanged=0;
        fnamestodelete={'subj1','subj2'};
        ofnfp=fullfile(todofiles(ind).folder,todofiles(ind).name);
        matcontents=load(ofnfp);
        for indz=fnamestodelete
            if isfield(matcontents,indz)
                matcontents=rmfield(matcontents,indz);
                ischanged=1;
            end
        end
        if ischanged
            save(ofnfp,'-struct','matcontents')
        end
        clear ofpfn matcontents
        
        continue
        
    elseif any(regexpi(todofiles(ind).name,'^README$|^.*.(m|txt)$'))              % m-files and TEXT files => search contents
        ofnfp=fullfile(todofiles(ind).folder,todofiles(ind).name);
        mcontents=fileread(ofnfp);

        updcontents=regexprep(mcontents,todoLUT(:,1),todoLUT(:,2));
        if ~strcmp(mcontents,updcontents)
            updcontents(updcontents == char(13)) = ''; % remove CR
    %         updcontents(updcontents == char(10)) = ''; % remove LF
            fid = fopen(ofnfp, 'wt+');
            fwrite(fid, updcontents);
            fclose(fid);
        end
        clear ofnfp mcontents updcontents
        continue

    else
        fprintf(1,'Unrecognized file "%s"\n',todofiles(ind).name)
    end
end
close(h)

for indz=1:size(todoLUT,1)
    token{1,indz}=randperm(16); % individual value for each folder/measurement
end

% SECOND RUN, rename filenames
h = waitbar(0,'Processing');
for ind=1:size(todofiles,1) 
    waitbar(ind/size(todofiles,1))   
    
      % run through all filenames
        for ind2=1:size(todoLUT,1)
             if any(regexp(todofiles(ind).name,todoLUT{ind2,1}))
                ofnfp=fullfile(todofiles(ind).folder,todofiles(ind).name); % old filename full path
                nfn=regexprep(todofiles(ind).name,todoLUT{ind2,1},todoLUT{ind2,2}); % new filename 
                if any(regexpi(todofiles(ind).name,'[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{2}'))    % Syngo XA61 data - StudyID hash
                    htemp=regexpi(todofiles(ind).name,'[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{2}','match');
                    htempn=htemp{1}(token{1});
                    nfn=regexprep(nfn,htemp,htempn); % new filename 
                end
                movefile(ofnfp,fullfile(todofiles(ind).folder,nfn))
                clear nfn ofpfn
                break % if found skip to next file
             end
        end
end
close(h)

% THIRD run, skip files, rename folders only
tododirs=dir(fullfile(pth,'**/*.*'));
tododirs(~[tododirs.isdir]) = []; 
tododirs=tododirs(~ismember({tododirs.name},{'.','..'}));
h = waitbar(0,'Processing');
for ind=1:size(tododirs,1)
    waitbar(ind/size(tododirs,1))
        for ind2=1:size(todoLUT,1)
            if any(regexp(tododirs(ind).name,todoLUT{ind2,1}))
                odnfp=fullfile(tododirs(ind).folder,tododirs(ind).name); % old dirname fullpath
                ndn=regexprep(tododirs(ind).name,todoLUT{ind2,1},todoLUT{ind2,2}); % new dirname
                movefile(odnfp,fullfile(tododirs(ind).folder,ndn))
                clear ndn odnfp
                break % if found skip to next dir
            end
        end

end
close(h)
