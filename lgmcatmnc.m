    datenow=datestr(now,'yyyymmdd');
   
        % create list of file prefixes
        timeser_files=dir(['*.0021600000.glob.nc']);
        
        for i=1:length(timeser_files)
            outfnm=[timeser_files(i).name(1:end-19),'.',datenow,'.glob.nc'];
            gluemnc_cat_timeseries(timeser_files(i).name(1:end-19),outfnm,0)
            
            % close any open files
            clear functions
            fclose('all');
        end
        
        % tidy up leftover timeseries fragments
        eval(['!find . -name "*.??????????.glob.nc" | xargs rm'])

    
