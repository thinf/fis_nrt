
dpath{1} = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\FSW1\data from deployment/';
dpath{2} = 'C:\Dropbox\Osci\FISP\#DATA\inductive system\FSE2\data from deployment/';

fn{1} ='microcat11_fulldata.txt';
fn{2} ='microcat12_fulldata.txt';
fn{3} ='microcat13_fulldata.txt';
fn{4} ='microcat14_fulldata.txt';
fn{5} ='microcat15_fulldata.txt';
fn{6} ='microcat16_fulldata.txt';

nme{1} = 'fsw1';
nme{2} = 'fse2';

for i = 1:numel(dpath)
    clear num t c p
    for j = 1:numel(fn);
        fname = [dpath{i} fn{j}];

        
%   22.0215,  0.00003,   -0.397, 09 Dec 2015, 14:00:01
%   22.5762,  0.00003,   -0.388, 09 Dec 2015, 16:00:01
%   22.6285,  0.00003,   -0.372, 09 Dec 2015, 18:00:01
%   22.8623,  0.00003,   -0.377, 09 Dec 2015, 20:00:01
%   22.8592,  0.00003,   -0.376, 09 Dec 2015, 22:00:01
%   22.6984,  0.00004,   -0.386, 10 Dec 2015, 00:00:01
%   22.4693,  0.00004,   -0.375, 10 Dec 2015, 02:00:01
   
         [t{j} c{j} p{j}, dd, mnd, yy, hh, mm, ss] =...
             textread(fname,'%f%f%f%s%s%s%s%s%s%*[^\n]',...
     'delimiter',',: ');
    % 'dd-mmm-yyyy HH:MM:SS'
    num{j} = nan([numel(dd),1]);
    for n = 1:numel(dd);
        num{j}(n) = datenum([dd{n} '-' mnd{n} '-' yy{n} ' ' hh{n} ':' mm{n} ':' ss{n}], 'dd-mmm-yyyy HH:MM:SS');
    end
    end
    cell2col(num,t,c,p)
    col2mat(num,t,c,p)
    
    make(nme{i}, 'num','t','c','p');
end

