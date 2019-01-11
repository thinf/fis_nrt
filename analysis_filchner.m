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
machine = 'remote_hh';
txtwrt = 0;
nrtfun = 0;
switch machine
    case 'local'
        dpath = 'C:\Dropbox\Osci\FISP\#DATA\SBD/';
        prfx = {'','','',''};
        workpath = 'C:\THINF\#TEMP\fis_nrt_tmp/';
        savepath = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\#IRIDIUM\processed/';
        have_jlab = 1;
    
    case 'remote'
        dpath = '/hs/datex/ingest/mooring_onisi/';
        prfx = {'fsw1/','fse2/','fne1/','fne2/'};
        workpath = '/home/csys/thatterm/fis_nrt/'; 
        savepath = '/home/csys/thatterm/Dropbox/AWIsync/SBD/';
        addpath('/home/csys/thatterm/mlab_lib/seawater');
        addpath('/home/csys/thatterm/mlab_lib/jlab');
        addpath('/home/csys/thatterm/mlab_lib/thinftools');
        addpath('/home/csys/thatterm/mlab_lib/eos/matlab');
        jlab_addpath
        have_jlab = 1;
        nrtfun = 1;
    case 'remote_hh'
        dpath = '/hs/datex/ingest/mooring/';
        prfx = {'fsw1/','fse2/','fne1/','fne2/'};
        workpath = '/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/fis_nrt/'; 
        savepath = '/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/data/';
        addpath('/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/mlab_lib/seawater');
        addpath('/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/mlab_lib/jlab');
        addpath('/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/mlab_lib/thinftools');
        addpath('/csys/ocean2/hhellmer/disk1/user1/work/projects/FISP/Mooring/mlab_lib/eos/matlab');
        jlab_addpath
        have_jlab = 1;
        nrtfun = 1;
end

stns = {'fsw1','fse2', 'fne1', 'fne2'};
imeis = {'/300234061031800_','/300234061032780_','/300234061035790_','/300234061030810_'};
nMCs = [6 6 5 5];
nADs = [4 4 3 3];
strtdates = [datenum(2015,12,24) datenum(2015,12,24) datenum(2016,12,1) datenum(2016,12,1)];
   
for ipi = 1:numel(stns)
    clearvars('-except','dpath','ipi', 'workpath', 'savepath', ...
        'machine', 'stns','imeis', 'nMCs', 'nADs','strtdates', 'prfx', ...
        'have_jlab','nrtfun','txtwrt')
    
    stn = stns{ipi}
    path = [dpath prfx{ipi} imeis{ipi}];
    
    nMC = nMCs(ipi); % number of Microcats
    nAD = nADs(ipi); % number of Aquadopps
    
    strtdate = strtdates(ipi);
    
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
        
        if txtwrt; dailySBD_filchner_txtwrite; end
        
        DayNo = DayNo + 1;
    end
    
    if ~isempty(broken_message)
    for i=1:numel(broken_message)
        delete(broken_message{i})
    end
     end
    
    %######## END OF DATAREADING SECTION, stacking and Saving ###############
    %%
    if have_jlab
        analysis_filchner_stack_datasave
        %  save([outpath stn '_' datestr(now,'yyyy_mm')],stn)
        analysis_filchner_stack_dataplot
    end
    if nrtfun
        analysis_filchner_nrt_plot
    end

end
