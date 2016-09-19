    %% stack Microcat data
    disp('stack and save microcat data')
    clear num p t c s th
    
    for n = 1:numel(data(1).Microcats)
        num{n} = [];
        p{n} = [];
        t{n} = [];
        c{n} = [];
        
        for i = 1:numel(data)
            hold on
            num_ = cell2mat(data(i).Microcats(n).Timestamp);
            
            p_ = data(i).Microcats(n).Pressure;
            t_ = data(i).Microcats(n).Temperature;
            c_ = data(i).Microcats(n).Conductivity;
            
            %abs(datenum(2016,1,1)-num)<500
            ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & p_ > 0 & p_ < 10000);
            
            num{n}=[num{n}; num_(ii)];
            p{n}=[p{n}; p_(ii)];
            t{n}=[t{n}; t_(ii)];
            c{n}=[c{n}; c_(ii)];
            
            if NumberOfDailyMessages==5
                num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
            end
        end
        
        
        %make dum p t th c s
        %mc{n} = dum;
    end
    cell2col(num,p,t,c);
    col2mat(num,p,t,c);
    s = sw_salt(10*c/sw_c3515,t,p);
    th = sw_ptmp(s,t,p,0);
    
    make mc num p t th c s
    
    %% stack Aquadopp data
    disp('stack and save Aquadopp data')
    clear num num_logger p t cv w head pitch roll amp1 amp2 amp3
    
    for n = 1:numel(data(1).Aquadopps)
        num{n} = [];
        num_logger{n} =[];
        p{n} = [];
        t{n} = [];
        cv{n} = [];
        w{n} = [];
        head{n} = [];
        pitch{n} = [];
        roll{n} = [];
        amp1{n} = [];
        amp2{n} = [];
        amp3{n} = [];
        
        for i = 1:numel(data)
            hold on
            %             for j = 1:numel(data(i).Aquadopps(n).Timestamp)
            %                 if isempty(data(i).Aquadopps(n).Timestamp{j})
            %                     data(i).Aquadopps(n).Timestamp{j} = nan;
            %                 end
            %             end
            
            %% latest version time:
            %if isempty( data(i).Aquadopps(n).Timestamp{1})
            %    data(i).Aquadopps(n).Timestamp{1} =nan;
            %end
            %num_ = data(i).Aquadopps(n).Timestamp{1}+[0:2:22]'/24; %cell2mat(data(i).Aquadopps(n).Timestamp);
            
            num_ = data(i).Aquadopps(n).Timestamp;
            num_ = [num_{:}]';
            
            num_logger_ = data(i).Aquadopps(n).LoggerTime;
            num_logger_ = [num_logger_{:}]';
            
%             %% scratch
%             n = 2;
%              num_ts =[]; 
%              num_lt = [];
%              for i = 1:numel(data);
%                  num_ = data(i).Aquadopps(n).Timestamp;
%                  num_ = [num_{:}]';
%                  num_ts = [num_ts; num_];
%                  
%                  num_ = data(i).Aquadopps(n).LoggerTime;
%                  num_ = [num_{:}]';
%                  num_lt = [num_lt; num_];
%              end
            %%
            
            
            p_ = data(i).Aquadopps(n).P;
            t_ = data(i).Aquadopps(n).T;
            cv_ = data(i).Aquadopps(n).U + 1i*data(i).Aquadopps(n).V;
            w_ = data(i).Aquadopps(n).W;
            head_ = data(i).Aquadopps(n).Head;
            pitch_ = data(i).Aquadopps(n).Pitch;
            roll_ = data(i).Aquadopps(n).Roll;
            amp1_ = data(i).Aquadopps(n).Amp1;
            amp2_ = data(i).Aquadopps(n).Amp2;
            amp3_ = data(i).Aquadopps(n).Amp3;
            
            %abs(datenum(2016,1,1)-num)<500
            %ii = find(isfinite((num_)) & abs(datenum(2016,1,1)-num_)<500 & abs(cv_)<1);
            try
            ii = find(isfinite((num_)) & abs(cv_)<1);
            dni = min(diff(num_(ii)));
            
            if ~isempty(ii) && numel(num{n})>0
            if num_(ii)< num{n}(numel(num{n}))
                
                disp(n)
                disp('found inconsistent time vector -- edges')
            elseif dni<0
                disp(n)
                disp('found inconsistent time vector -- inside')
                
                
            end
            end
             catch
                keyboard
            end
           
            num_logger{n}=[num_logger{n}; num_logger_(ii)];
            num{n}=[num{n}; num_(ii)];
            p{n}=[p{n}; p_(ii)];
            t{n}=[t{n}; t_(ii)];
            cv{n}=[cv{n}; cv_(ii)];
            w{n}=[w{n}; w_(ii)];
            head{n}=[head{n}; head_(ii)];
            pitch{n}=[pitch{n}; pitch_(ii)];
            roll{n}=[roll{n}; roll_(ii)];
            amp1{n}=[amp1{n}; amp1_(ii)];
            amp2{n}=[amp2{n}; amp2_(ii)];
            amp3{n}=[amp3{n}; amp3_(ii)];
            
            if NumberOfDailyMessages==5
                num{n}(num{n}>datenum(2016,1,21))=datenum(2016,1,21);
            end
        end
        
        
        %make dum p t th c s
        %mc{n} = dum;
    end
    cell2col(num, num_logger,p,t,cv,w,head,pitch,roll,amp1,amp2,amp3);
    col2mat(num, num_logger,p,t,cv,w,head,pitch,roll,amp1,amp2,amp3);
  
    make ad num num_logger p t cv w head pitch roll amp1 amp2 amp3
    
%     %% stack Housekeeping data
%     disp('stack and save Housekeeping data')
%     clear num p t cv w head pitch roll amp1 amp2 amp3
%     
%     for n = 1:numel(data(1).Aquadopps)
%     end
    
    timestamp = now;
    make dum mc ad timestamp files path
    eval([stn '= dum;'])
    save([outpath stn '_' datestr(now,'yyyy_mm')],stn)

    %######## END OF DATASAVING SECTION, make standard plots ###############
