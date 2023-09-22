%{
 This script takes as input individual monthly files containing B-grid velocities:
    U-vel (with dimensions time (length 1), st_ocean, yu_ocean, xu_ocean),
    V-vel (with dimensions time (length 1), st_ocean, yu_ocean, xu_ocean),
    W-vel (with dimensions time (length 1), sw_ocean, yt_ocean, xt_ocean)
    
    and resaves the files in a different directory into monthly files containing u,v,w on
    an A-grid with dimensions:
    time (length 1), st_ocean (+1 length of above files for bottom boundary condition), yu_ocean, xu_ocean 
    
    WW: Also replaces all land masks by zero velocities. 
%}

exp = 'hybrid_v5_rel04_BC5_ne120_t12_pop62'
startyear = 66;
startmonth = 1;

% Path to velocity data
path2files = strcat('/usr/projects/cesm/ASD/for_CMS_SO/');

% Path to saving directory
path2save = strcat('/usr/projects/cesm/ASD/for_CMS_SO_Edited/')

files = dir(path2files);
files = sort({files.name});

% Go through each different file ending with 'u.nc' to get the date prefix
% and then extract u,v,w from that 5-day period. Make sure to extract the
% u-cell mask.

for q=1:length(files)
    file = char(files(q));
    if length(file) < 5
        continue
    end
    if ~(strcmp(file(end-3:end),'u.nc'))
        continue
    end
    year = str2num(file(8:11))
    if year > 66
        continue
    end
    month = str2num(file(12:13));
    day = str2num(file(14:15));
    disp(['Year: ',num2str(year),' Month: ',num2str(month),' Day: ',num2str(day)])
    u_file = strcat(path2files,file);
    v_file = strcat(path2files,file(1:end-4),'v.nc');
    w_file = strcat(path2files,file(1:end-4),'w.nc');
    % Take dimension values and masks only once to cut I/O time.
    if ((year == startyear) && (month == startmonth))

        time = ncread(u_file,'time');
%WW	convert from cm to m
        st_ocean = 0.01*ncread(u_file,'z_t');
        sw_ocean = 0.01*ncread(u_file,'z_w');
%WW	extract one longitude/latitude line from the 2D ULAT and ULONG fields
        dum = ncread(u_file,'ULAT');yu_ocean = dum(3000,:);
        dum = ncread(u_file,'ULONG');xu_ocean = dum(:,500);
%WW	Make some adjustments to longitude to make it monotonic
        xx = find(xu_ocean(2001:end) < 0);xu_ocean(2000+xx) = xu_ocean(2000+xx)+360;
%WW	extract one longitude/latitude line from the 2D TLAT and TLONG fields
        dum = ncread(u_file,'TLAT');yt_ocean = dum(3000,:);
        dum = ncread(u_file,'TLONG');xt_ocean = dum(:,500);
%WW	Make some adjustments to longitude to make it monotonic
        xx = find(xt_ocean(1:2000) > 180);xt_ocean(xx) = xt_ocean(xx)-360;

        % Need buffer latitude layers above 30°S for 1) Interpolation of W and 2) For
        % the way CMS works.
        lat_n_index = length(yu_ocean)-1;
%       lat_n_index = lat_n_index + 2;

        umask = ncread(u_file,'UVEL'); umask = umask(:,1:lat_n_index,:);
        vmask = ncread(v_file,'VVEL'); vmask = vmask(:,1:lat_n_index,:);

%WW  	fiddle with the masks a little. Masks are either indicated by NaN, -1, or possibly zero, if masking has not been applied
%WW	consistently.
        umask(isnan(umask))=1;
        umask(umask==-1)=1;
        umask(umask==0)=1;
        umask(umask~=1)=0;

        vmask(isnan(vmask))=1;
        vmask(vmask==-1)=1;
        vmask(vmask==0)=1;
        vmask(vmask~=1)=0;

%WW	A few inconsistencies here and there between UVEL and VVEL based masks.
%WW	Make sure that 'ocean' is the larger subset.
        umask(vmask==0)=0;
        
%WW	Same for w mask
        wmask = ncread(w_file,'WVEL'); wmask = wmask(:,1:lat_n_index+1,:);
        wmask(isnan(wmask))=1;
        wmask(wmask==-1)=1;
        wmask(wmask==0)=1;
        wmask(wmask~=1)=0;

        yu_ocean = yu_ocean(1:lat_n_index);
        yt_ocean = yt_ocean(1:lat_n_index);
        diff_bottom = st_ocean(end)-sw_ocean(end); % Need this for later
        st_ocean(end+1)=st_ocean(end)+diff_bottom*2.0;

    end
%WW	Rescale from cm/s to m/s
    u = 0.01*ncread(u_file,'UVEL'); u = u(:,1:lat_n_index,:);
    v = 0.01*ncread(v_file,'VVEL'); v = v(:,1:lat_n_index,:);
    w = 0.01*ncread(w_file,'WVEL'); w = w(:,1:lat_n_index+1,:);
    % Extra buffer north for w because we lose a layer from horizontal
    % interpolation
    
    n_xt = length(xt_ocean);
    n_yt = length(yt_ocean);
    n_sw = length(sw_ocean);
    
    % Add a bottom layer that will contain a boundary so that CMS has 0 velocities to interpolate with
    % Take away mask, add extra level (all masked), then reapply the mask.
%WW	I slightly changed this. I'm setting all the masked-out regions to zero
    u_null_mask = u;
    u_null_mask(umask == 1) = 0;
    u_null_mask = cat(3,u_null_mask,zeros(n_xt,n_yt));
    u_regridded=u_null_mask;
    
    % Add a bottom layer that will contain a boundary so that CMS has 0 velocities to interpolate with
    % Take away mask, add extra level (all masked), then reapply the mask.
    v_null_mask = v;
    v_null_mask(umask == 1) = 0;
    v_null_mask = cat(3,v_null_mask,zeros(n_xt,n_yt));
    v_regridded=v_null_mask;
    
    % Replace mask with large values so interpolation gives large values if one of the 4 points is land [(large + 0 + 0 + 0)/4 = still large is the idea here]

%WW	 Again, using zeros for the masked values.

    w_null_mask = w(:,:,:);
    w_null_mask(wmask == 1) = 0;

    % Horizontal Interpolation

%WW	The ordering is differentin matlab, so circshift has to be applied to the first dimension; and -1 will shift it to the east. 

    w_sw = w_null_mask(:,1:end-1,:);
    w_nw = w_null_mask(:,2:end,:);
    w_se = circshift(w_sw,[-1 0 0]); % Shift values to the east by one
    w_ne = circshift(w_nw,[-1 0 0]); % Shift values to the east by one
    w2d = 0.25*(w_sw+w_nw+w_se+w_ne);
    
    % Vertical Interpolation

%WW	For POP, w(1) is at the surface. It is the bottom value that is not explicitly included and is implictly assumed to be zero.

    w_regridded = zeros(n_xt,n_yt,n_sw);
    for k=1:n_sw-1
      dzk=st_ocean(k)-sw_ocean(k);
      dz=sw_ocean(k+1)-sw_ocean(k);
      w_regridded(:,:,k) = w2d(:,:,k) + (w2d(:,:,k+1)-w2d(:,:,k)).*dzk./dz;
    end
    k=n_sw;
    dzk=st_ocean(k)-sw_ocean(k);
    dz=2*dzk;
    w_regridded(:,:,k) = w2d(:,:,k) - w2d(:,:,k).*dzk./dz;

    % Insert layer at the bottom (because of trilinear interpolation method in CMS)
    w_regridded = cat(3,w_regridded,zeros(n_xt,n_yt));

%WW	Apply umask to the w field

    w_regridded(umask == 1) = 0;

    % Make the cell below bottom cells have negative the velocities of the cell above. Then, halfway between the two (approx the sw_ocean bottom of the domain) will have W-vel = 0.
    for i = 1:n_xt
        for j = 1:n_yt
            xx=find(umask(i,j,:)==1,1);
            if xx > 1
              w_regridded(i,j,xx) = -w_regridded(i,j,xx-1);
            end
        end
    end

    % Only w is regridded but the landmasks were changed to 0 for all three fields.
    uvar(:,:,:,1) = u_regridded;
    vvar(:,:,:,1) = v_regridded;
    wvar(:,:,:,1) = w_regridded;

    %WW: Now prepare and write netcdf file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save u,v,w on A-grid to same file for given 5-day period

    % Name the file
    outfile_uvw = strcat(path2save,file(1:end-4),'.nc');

    % Create writable netcdf file at location outfile_uvw
    ncFileID = netcdf.create(outfile_uvw,'NETCDF4');

    % Create dimensions
    % Record dimension of time with length 'None'
%WW	Didn't recognize 'none', so I made it equal to 1?
    timeDimID=netcdf.defDim(ncFileID,'time',1)

    depthDimID=netcdf.defDim(ncFileID,'st_ocean', length(st_ocean))
    latDimID=netcdf.defDim(ncFileID,'yu_ocean', length(yu_ocean))
    lonDimID=netcdf.defDim(ncFileID,'xu_ocean', length(xu_ocean))

    % Create variables
    ti = netcdf.defVar(ncFileID,'time', 'NC_FLOAT', timeDimID);
    st = netcdf.defVar(ncFileID,'st_ocean', 'NC_FLOAT', depthDimID);
    yu = netcdf.defVar(ncFileID,'yu_ocean', 'NC_FLOAT', latDimID);
    xu = netcdf.defVar(ncFileID,'xu_ocean', 'NC_FLOAT', lonDimID);

%   u_var = netcdf.defVar(ncFileID,'u', 'NC_FLOAT', [depthDimID latDimID lonDimID timeDimID]);
%   v_var = netcdf.defVar(ncFileID,'v', 'NC_FLOAT', [depthDimID latDimID lonDimID timeDimID]);
%   w_var = netcdf.defVar(ncFileID,'w', 'NC_FLOAT', [depthDimID latDimID lonDimID timeDimID]);
    u_var = netcdf.defVar(ncFileID,'u', 'NC_FLOAT', [lonDimID latDimID depthDimID timeDimID]);
    v_var = netcdf.defVar(ncFileID,'v', 'NC_FLOAT', [lonDimID latDimID depthDimID timeDimID]);
    w_var = netcdf.defVar(ncFileID,'w', 'NC_FLOAT', [lonDimID latDimID depthDimID timeDimID]);

%WW	Disable fill value, so that CMS will see only zeros.
    netcdf.defVarFill(ncFileID,u_var,true,0)
    netcdf.defVarFill(ncFileID,v_var,true,0)
    netcdf.defVarFill(ncFileID,w_var,true,0)

% Exit define mode
    netcdf.endDef(ncFileID);

% Write the variables to file
    netcdf.putVar(ncFileID,ti,time)
    netcdf.putVar(ncFileID,st,st_ocean)
    netcdf.putVar(ncFileID,yu,yu_ocean)
    netcdf.putVar(ncFileID,xu,xu_ocean)

    netcdf.putVar(ncFileID,u_var,uvar)
    netcdf.putVar(ncFileID,v_var,vvar)
    netcdf.putVar(ncFileID,w_var,wvar)

% Re-enter define mode
    netcdf.reDef(ncFileID);

% Add the attributes
    netcdf.putAtt(ncFileID,ti,'units','days since 0001-01-01 00:00:00');
    netcdf.putAtt(ncFileID,ti,'calender_type', 'JULIAN');
    netcdf.putAtt(ncFileID,ti,'calender','JULIAN');

    netcdf.putAtt(ncFileID,st,'units','meters');
    netcdf.putAtt(ncFileID,st,'long_name','tcell depth');
    netcdf.putAtt(ncFileID,st,'positive','down');

    netcdf.putAtt(ncFileID,yu,'units','degrees N');
    netcdf.putAtt(ncFileID,yu,'long_name','ucell latitude');

    netcdf.putAtt(ncFileID,yu,'units','degrees E');
    netcdf.putAtt(ncFileID,yu,'long_name','ucell longitude');

    netcdf.putAtt(ncFileID,u_var,'units','m/s');
    netcdf.putAtt(ncFileID,yu,'long_name','zonal velocity');

    netcdf.putAtt(ncFileID,u_var,'units','m/s');
    netcdf.putAtt(ncFileID,yu,'long_name','meridional velocity');

    netcdf.putAtt(ncFileID,u_var,'units','m/s');
    netcdf.putAtt(ncFileID,yu,'long_name','vertical velocity');

    netcdf.close(ncFileID);

end


        
