% First attemtp to read and plot the FISP Iridium data
% still heavily under construction
% thinf, 04.02.2015, tore.hattermann@awi.de

clear all;
%close all;

% Change this path to the folder and filename prefix that the SBD messages
% are stored in.
% Please note that datalogger at FSW1 SBD messages begin with 300234061031800_
% Please note that datalogger at FSE2 SBD messages begin with 300234061038780_
%
% path = 'E:\SBDs\from Gerd\300234061031800_';
% path = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\800/300234061031800_';
paths{1} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061031800_';  stns{1} = 'FSW1';
paths{2} = 'C:\Dropbox\Osci\FISP\#DATA\SBD/300234061032780_'; stns{2} = 'FSE2';
for ipi =  1:numel(paths)
    clearvars('-except','paths','ipi', 'stns')
    path = paths{ipi}
    stn = stns{ipi}
    % specify folder where temporary files will be created. Also the files
    % "blankmsg1.sbd" to "blankmsg7.sbd" MUST be stored here.
    outpath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\processed/';
    workpath = 'C:\THINF\#TEMP\fis_nrt_tmp';
    
    nMC = 6; % number of Microcats
    nAD = 4; % number of Aquadopps
    
    % SBD file extension
    filetype = '.sbd';
    % location of blankmessages
    blankmsgs = ['./', 'blankmsg'];
    
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
        
        [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
            = dailySBD_filchner(nMC, nAD, fids, workpath, msg_format);
        
        %[data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
        %    = dailySBD_filchner(nMC, nAD, fids, outpath, add_dummybytes);
        
        DayNo = DayNo + 1;
    end
    %######## END OF DATAREADING SECTION, stacking and Saving ###############
    %%
 analysis_filchner_stack_datasave
 %  save([outpath stn '_' datestr(now,'yyyy_mm')],stn)
 analysis_filchner_stack_dataplot

end
