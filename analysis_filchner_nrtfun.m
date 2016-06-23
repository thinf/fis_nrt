function analysis_filchner_nrtfun(msg_prefix,stn,workpath,savepath)
% First attemtp to read and plot the FISP Iridium data
% still heavily under construction
% thinf, 04.02.2015, tore.hattermann@awi.de
% thinf, 14.04.2016. extended to run as standalone

% Janik:
% hier die beiden Pfäde:
%  
% /hs/datex/ingest/mooring/fispocean/fsw1/data/
% /hs/datex/ingest/mooring/fispocean/fse2/data/
%  
% Es sind nun 2 XML Dateien dafür zuständig. 
% 
% Das Skript liegt aktuell /hs/datex/ingest/mooring/fispocean/ dort. 
% 
%
%
% Change this path to the folder and filename prefix dthat the SBD messages
% are stored in.
% Please note that datalogger at FSW1 SBD messages begin with 300234061031800_
% Please note that datalogger at FSE2 SBD messages begin with 300234061038780_
%
% path = 'E:\SBDs\from Gerd\300234061031800_';
% path = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\800/300234061031800_';
%paths{1} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061031800_';  stns{1} = 'FSW1'
%paths{2} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061032780_'; stns{2} = 'FSE2';

%paths{1} = 'data/300234061031800_';  stns{1} = 'FSW1'
%paths{2} = 'data/300234061032780_'; stns{2} = 'FSE2';

% for srcipt only
%for ipi = 1:numel(paths)
%    clearvars('-except','paths','ipi', 'stns')
%    path = paths{ipi}
%    stn = stns{ipi}

% specify folder where temporary files will be created. Also the files
    % "blankmsg1.sbd" to "blankmsg7.sbd" MUST be stored here.
    % workpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\processed/';
    % workpath = './';
    
    % workpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\nrt_demo\';
    %  workpath = '/hs/datex/ingest/mooring/fispocean/';
    
    nMC = 6; % number of Microcats
    nAD = 4; % number of Aquadopps
    
    
    % SBD file extension
    path = msg_prefix; %search path for messages
    filetype = '.sbd';
    blankmsgs = ['./', 'blankmsg.sbd'];
    
    NumberOfDailyMessages = 7;
    % The Filchner dataloggers have 7 SBD messages sent by Iridium every day.
    % Earlier SBDs (ca. before 21.01.2016) from IMEI ...800 had 5 messages only.
    % The program can read these messages as well, but the Aquadopp data from this
    % period will make no sense and the Microcat data will have gaps.
    % To properly read those older files, you can set the value back to 5, but
    % then the Program will fail to read newer files
    
    if NumberOfDailyMessages == 7;
        msg_format = 1; % standard Filchner south data format
    elseif NumberOfDailyMessages == 5;
        msg_format = 0; % set to 0 for old AD-format/
    else
        error('STOP, number of daily messages suggest unknown dataformat')
    end
    
    %% ######## END OF USER DEFINED SECTION ###############
      
    analysis_filchner_SBDsetup
    
    %Loop through messages and process data.
    DayNo = 1;
    while MessageNo < (MaxNumberOfDays * NumberOfDailyMessages) && DayNo <= MaxNumberOfDays;
        
        analysis_filchner_SBDloop
        
        %Strip data from files and put in 'data' structure.
        %     [data(DayNo).Housekeeping, data(DayNo).Microcat1, data(DayNo).Microcat2,...
        %         data(DayNo).Microcat3, data(DayNo).Microcat4, data(DayNo).Microcat5,...
        %         data(DayNo).Microcat6, data(DayNo).Aquadopp1, data(DayNo).Aquadopp2,...
        %         data(DayNo).Aquadopp3, data(DayNo).Aquadopp4]...
        %         = dailySBD_filchner_FSW1_th(MessageArray{1, 3, DayNo},...
        %         MessageArray{2, 3, DayNo}, MessageArray{3, 3, DayNo},...
        %         MessageArray{4, 3, DayNo}, MessageArray{5, 3, DayNo}, workpath);
        %
        
        % [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
        %         = dailySBD_filchner_FSW1_th(nMC, nAD, MessageArray{1, 3, DayNo},...
        %         MessageArray{2, 3, DayNo}, MessageArray{3, 3, DayNo},...
        %         MessageArray{4, 3, DayNo}, MessageArray{5, 3, DayNo}, workpath);
        
        [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
            = dailySBD_filchner(nMC, nAD, fids, workpath, msg_format);
        
        % write out txt data for NRT database
        % uses "stn" and "workpath" from input
        dailySBD_filchner_txtwrite
        
        %%
        DayNo = DayNo + 1;
    end
% end % run over paths for script only
end
