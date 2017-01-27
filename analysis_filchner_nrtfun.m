function analysis_filchner_nrtfun(msg_prefix,stn,workpath,savepath)

path = msg_prefix;
    
    if strcmp(stn,'fsw1') || strcmp(stn,'fse2')
        strtdate = datenum(2015,12,24);
        nMC = 6; % number of Microcats
        nAD = 4; % number of Aquadopps

    elseif strcmp(stn,'fne1') || strcmp(stn,'fne2')
        strtdate = datenum(2016,12,1);
        nMC = 5; % number of Microcats
        nAD = 3; % number of Aquadopps

    end 
    
    
    % specify folder where temporary files will be created. Also the files
    % "blankmsg1.sbd" to "blankmsg7.sbd" MUST be stored here.
    %
    % workpath = './';
    
    % specify folder where figures / processed data are saved
    %
    % outpath = './';
    
    
    
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
   
    % analysis_filchner_SBDsetup
     % %MaxNumberOfDays - this number should be changed to be equal to the number
    % %of days data you want to process.  It might be worth changing it to (say)
    % %365 days now, although this will mean the script takes a long time to
    % %run.
    % now determined automatically:
    MaxNumberOfDays = round(now-strtdate);%100;
    files = dir([path '*']);
    
    %set up 'message array' to store date, filename and fileID
    MessageArray = cell(NumberOfDailyMessages, 3, MaxNumberOfDays);
    
    % initialize output variable that will contain the actual data
    data = struct([]);
    broken_message = struct([]);
    bm = 0;
    
    %% loop to find first valid file in directory.  This means that any files
    %that were sent by the modem and then deleted before deployment (eg test
    %messages etc) are ignored.
    for MessageNo=0:(MaxNumberOfDays*NumberOfDailyMessages)
        
        testFilename = strcat(path, sprintf('%06d', MessageNo),filetype);
        
        if fopen(testFilename) ~= -1;
            disp('start with message')
            disp(testFilename)
            break
        end
        
    end
    
    %Loop through messages and process data.
    DayNo = 1;
    while MessageNo < (MaxNumberOfDays * NumberOfDailyMessages) && DayNo <= MaxNumberOfDays;
        
        analysis_filchner_SBDloop
        
        [data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
            = dailySBD_filchner(nMC, nAD, fids, workpath, msg_format);
        
        %[data(DayNo).Housekeeping, data(DayNo).Microcats,data(DayNo).Aquadopps]...
        %    = dailySBD_filchner(nMC, nAD, fids, outpath, add_dummybytes);
        
       dailySBD_filchner_txtwrite
        
        DayNo = DayNo + 1;
    end
    
    if ~isempty(broken_message)
    for i=1:numel(broken_message)
        delete(broken_message{i})
    end
     end
    
    %######## END OF DATAREADING SECTION, stacking and Saving ###############
   exit
end