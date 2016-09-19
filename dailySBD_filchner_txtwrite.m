%% write out data of the day
% producing txtfiles suitable for the
% AWI near real-time database
if isfinite(data(DayNo).Housekeeping.Date)
    % dfile = [paths{ipi} stns{ipi} '_' datestr(datenum(data(DayNo).Housekeeping.Date),'yyyy-mm-dd')];
    dfile = [savepath stn '_' datestr(datenum(data(DayNo).Housekeeping.Date),'yyyy-mm-dd')];
    %%
    if ~exist([dfile '.hdr'],'file')
        diary([dfile '.hdr'])
        disp(data(DayNo).Housekeeping)
        diary off
        
        
        %%  Write Microcat data
        fields = {'Timestamp','Pressure','Temperature','Conductivity'};
        units = {'JulianDays','dBar','degC','mS/cm'};
        suffix = '_Microcats';
        file = [dfile suffix '.txt'];
        fid = fopen(file,'w+');
        dmat = [];
        for mcs = 1:nMC;
            % write header
            
            for iif = 1:numel(fields)
                fprintf(fid,'%s',['MC' num2str(mcs) ':' fields{iif} ':' units{iif} ';']);
            end
            dmat = [dmat, [data(DayNo).Microcats(mcs).Timestamp{:}]',...
                data(DayNo).Microcats(mcs).Pressure, ...
                data(DayNo).Microcats(mcs).Temperature, ...
                data(DayNo).Microcats(mcs).Conductivity*10,...
                ];
            
        end
        
        fprintf(fid,'\n');
        fclose(fid);
        
        dlmwrite(file,dmat,'delimiter',';','precision',9,'-append')
        
        %%  Write Aquadopp data
        
        
        fields = {'Timestamp','U','V','W','Head','Pitch','Roll','Pressure','Temperature'};
        units = {'JulianDays','m/s','m/s','m/s','deg','deg','deg','dBar','degC'};
        suffix = '_Aquadopps';
        file = [dfile suffix '.txt'];
        fid = fopen(file,'w+');
        dmat = [];
        
        for mcs = 1:nAD;
            % write header
            
            for iif = 1:numel(fields)
                fprintf(fid,'%s',['AD' num2str(mcs) ':' fields{iif} ':' units{iif} ';']);
            end
            
            adtime = [data(DayNo).Aquadopps(mcs).Timestamp{:}]';
            
            %if isempty( data(DayNo).Aquadopps(mcs).Timestamp{1})
            %   data(DayNo).Aquadopps(mcs).Timestamp =-9999;
            %end
            %num_ = data(i).Aquadopps(n).Timestamp{1}+[0:2:22]'/24; %cell2mat(data(i).Aquadopps(n).Timestamp);
            adtime(isnan(adtime))=-9999;
            
            dmat = [dmat, adtime,...
                data(DayNo).Aquadopps(mcs).U, ...
                data(DayNo).Aquadopps(mcs).V, ...
                data(DayNo).Aquadopps(mcs).W, ...
                data(DayNo).Aquadopps(mcs).Head, ...
                data(DayNo).Aquadopps(mcs).Pitch, ...
                data(DayNo).Aquadopps(mcs).Roll, ...
                data(DayNo).Aquadopps(mcs).P, ...
                data(DayNo).Aquadopps(mcs).T, ...
                ];
        end
        
        fprintf(fid,'\n');
        fclose(fid);
        dlmwrite(file,dmat,'delimiter',';','precision',9,'-append')
    end
end