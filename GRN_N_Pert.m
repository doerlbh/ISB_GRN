% GRN_N_Pert
% Perturbating n-gene multistate networks
% Author: Baihan Lin
% Date:   August 2016


clear all;
close all;

rng(111);                 % randomizer

pathN = '/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/';
pathI = '/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/';
system(['mkdir ' pathN]);

% prompt = 'What is your folder?: ';
% disp('e.g.  /Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');
% pathI = input(prompt,'s');

% [~,list] = system(['find ' path ' -type f -name "MX*.txt"']);
%
% files = strsplit(list);
% length(files);
%
% for countfile = 1:length(files)
%     trials{countfile} = files{countfile}(1:end-11);
% end
%
% trials = unique(trials);
% trials = trials(~cellfun('isempty',trials));
%
% index = strfind(files, trials{1});
%
% for t = 1 : length(trials)
%     filename = trials{t};
%     indexsep = strfind(filename,'/');
%     last = indexsep(end);
%
%     index = strfind(files, filename);
%     indexmat = cell2mat(index);
%     [~,first,~] = unique(indexmat, 'first');
%     fast_arrange_Intan_RHD(files{first(1)});
%
%     casename = strtok(filename(last+1:end), '_');
%

P = dlmread('/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');
M = dlmread('/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');
A = dlmread('/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');
N = dlmread('/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');
Nss = dlmread('/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest11-20160818/');

ph = 100;
y0 = Nss(t,:);

Pt = P;
Nt = N;
At = A;
Mt = M;
P(3,4) = 100*P(3,4);

y = pertrun(pet,0,n,P,N,A,M,ph,1,y0,sstrshd,zthrs,Pt,Nt,At,Mt);
    
parsave(strcat(pathN, fname, '-y.txt'), y, '-ascii');



