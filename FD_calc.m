%% CALCULATION OF FRAMEWISE DISPLACEMENT

function FD=FD_calc(mov_par,fd_type)
    if ~exist("fd_type","var")
        fd_type = 1;  
    end
    
    %Regressors
    reg_trans     = mov_par(:,1:3);
    reg_rot       = mov_par(:,4:6);
    dif_trans     = [0 0 0; (reg_trans(2:end,:)-reg_trans(1:(end-1),:))];
    dif_rot       = [0 0 0; (reg_rot(2:end,:)-reg_rot(1:(end-1),:))];


    %Framewise Displacement
    switch fd_type
        case 1 % Jonathan Power
            abs_trans = abs(dif_trans);
            abs_rot   = abs(dif_rot);
            % Rotational displacements were converted from degrees to
            % millimeters by calculating displacement on the surface of a
            % sphere of radius 50 mm, which is approximately the mean distance
            % from the cerebral cortex to the center of the head. - z èlánku
            % Power 2012 Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion
            fd = abs_trans(:,1)+abs_trans(:,2)+abs_trans(:,3)+50*pi/180*(abs_rot(:,1)+abs_rot(:,2)+abs_rot(:,3)); %Degrees are here transformed to radians
        case 2 % Chao-Gan Yan
            fd = sqrt((dif_trans(:,1).^2)+(dif_trans(:,2).^2)+(dif_trans(:,3).^2));
        case 3 % Koene Van Dijk
            fd = sqrt((reg_trans(2:end,1).^2)+(reg_trans(2:end,2).^2)+(reg_trans(2:end,3).^2))-sqrt((reg_trans(1:(end-1),1).^2)+(reg_trans(1:(end-1),2).^2)+(reg_trans(1:(end-1),3).^2));
            fd = [0; fd];
    end
FD=fd;




