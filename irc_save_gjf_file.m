function irc_save_gjf_file(ms0, odir, gjf_odir, fullgtemplname)
    order=1:ms0.atomnum;
    if exist([odir filesep gjf_odir filesep ms0.desc '.gjf'],'file')~=2 %create input Gaussian file
        if exist([odir filesep gjf_odir],'dir')~=7
            mkdir(odir,gjf_odir);
        end
        savemolgs([odir filesep gjf_odir],ms0,3,order,fullgtemplname); %Gaussian with XYZ
    end
end