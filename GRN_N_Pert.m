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

    P = = dlmread(strcat(pathI,));
    M = = dlmread();
    A = Ar*diag(P)/Pr;
    N = zeros(n,n);
    
    for node = 1:n
        N(node,:) = randsample(-1:1,n,true);
        N(node,node) = 0;
    end
    
    %     xr = zeros(trial,n,1);
    
    Nss = [];
    %     t = 1;
    for t = 1:trial
        
        ph = 100;
        y0 = x0(t,:);
        [y,ph,ss] = netrun(n,P,N,A,M,ph,1,y0,tend,sstrshd,zthrs);
        if ss == 1
            ssn = y(size(y,1),:);
%             if ~prod((abs(ssn)<zthrs))
                
                if size(Nss,1) > 0
                    for st = 1:size(Nss,1)
                        samess(count) = prod((Nss(st,:)-ssn)<ddtrshd);
                    end
                else
                    samess(count) = 0;
                end
                
                if samess(count) == 0
                    Nss = vertcat(Nss, ssn);
                    
                    fig = figure;
                    plot(y);
                    title(strcat('X',num2str(count),'-N',num2str(n),'-T',num2str(t)));
                    filename = strcat(pathN, 'X',num2str(count),'-N',num2str(n),'-T',num2str(t),'.png');
                    parsaveas(gcf, filename,'png');
                    close gcf;   
                end
%             end
        end
    end
    ssc = size(Nss,1);
    
    if ssc > 1
        p = p + 1;
        fname = strcat('MX',num2str(count),'-N',num2str(n),'-C',num2str(ssc));
    else
        fname = strcat('UX',num2str(count),'-N',num2str(n),'-C',num2str(ssc));
    end
    
    parsave(strcat(pathN, fname, '-Nss.txt'), Nss, '-ascii');
    parsave(strcat(pathN, fname, '-P.txt'), P, '-ascii');
    parsave(strcat(pathN, fname, '-A.txt'), A, '-ascii');
    parsave(strcat(pathN, fname, '-N.txt'), N, '-ascii');
    parsave(strcat(pathN, fname, '-M.txt'), M, '-ascii');
    


