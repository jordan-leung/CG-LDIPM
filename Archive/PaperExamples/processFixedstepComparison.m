clc
clear all
close all

% In format: (n,m,condNum)
CASE_DATA = [100 200 100;...
    100 1000 100;
    500 1000 100;
    500 2000 100;
    1000 2000 100;
    1000 5000 100];
NCases = size(CASE_DATA,1);

%% Process the data

% Initialize
CG_ITER_MEAN = zeros(NCases,2);
CG_ITER_STD = zeros(NCases,2);
LDIPM_ITER_MEAN = zeros(NCases,2);
LDIPM_ITER_STD = zeros(NCases,2);
INFEAS_NUM = zeros(NCases,1);

% Load the data into a structure
for i = 1:NCases
    n = CASE_DATA(i,1);
    m = CASE_DATA(i,2);
    condNum = CASE_DATA(i,3);
    loadStr = ['./Data/fixedstepData_n',num2str(n),'_m',num2str(m),'_cond',num2str(condNum)];
    load(loadStr);
    if i == 1
        DATA = repmat(saveData,NCases,1);
    else
        DATA(i) = saveData;
    end
    CG_ITER_MEAN(i,:) = mean(saveData.CG_ITERS)/1e3;
    CG_ITER_STD(i,:) = std(saveData.CG_ITERS)/1e3;
    LDIPM_ITER_MEAN(i,:) = mean(saveData.LDIPM_ITERS);
    LDIPM_ITER_STD(i,:) = std(saveData.LDIPM_ITERS);
    INFEAS_NUM(i) = sum(~saveData.FEASFLAG);
end



percent = (CG_ITER_MEAN(:,1) - CG_ITER_MEAN(:,2))./CG_ITER_MEAN(:,1)*100
min(percent)
max(percent)

AVG = CG_ITER_MEAN*1e3./(LDIPM_ITER_MEAN)