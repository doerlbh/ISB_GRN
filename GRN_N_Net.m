% GRN_N_Net
% Creating n-gene multistate networks
% Author: Baihan Lin
% Date:   August 2016


clear all;
close all;

rng(100);                 % randomizer

pathN = '/Users/DoerLBH/Dropbox/git/ISB_GRN/data/Mtest15-20160818/';
system(['mkdir ' pathN]);

num = 20;         % network number
n = 9;           % gene
trial = 50;      % trial
ddtrshd = 1;       % threshold for different states
sstrshd = 1;     % threshold for equilibrium of steady states
tend = 3000;     % threshold for equilibrium vs. non-equlibrium
zthrs = 1;       % threshold for non-trivial states

% Starting values
Nr = 50;  % range of initial states
Pr = 10;  % range of parameters
Ar = 100; % range of synthesis rate
pow = 6;  % range of power increase

x0 = Nr*rand(trial, n);

samess = zeros(num,1);

% Direct parameter declaration:
x = zeros(n,1);   % gene number
A = zeros(n,n);   % gene sythesis rate
M = zeros(n,1);   % gene degradation rate
P = zeros(n,n);   % parameters, i.e. gene-gene interaction (repress or activate)

p = 1;           % num of generated networks (starting 1)

% count = 1;
parfor count = 1:num
    
    P = Pr*rand(n, n);
    M = rand(n,1);
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
        [y,ph,ss] = netrun(n,P,N,A,M,ph,1,y0,tend,sstrshd,zthrs,pow);
        if ss == 1
            ssn = real(y(size(y,1),:));
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
                    plot(real(y));
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
    
end

