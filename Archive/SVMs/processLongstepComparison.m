clc
clear all
close all


%% Process the data


clc
clear all
close all

% Specify the folder where the files live.
myFolder = '/Users/jordan/Documents/GitHub/CG-LDIPM/MPC_Examples/Randomized/Data/LongStep';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% Preassign table arrays
NCases = length(theFiles);

% Initialize
SIZES = zeros(NCases,4);
CG_ITER_MEAN = zeros(NCases,2);
CG_ITER_STD = zeros(NCases,2);
LDIPM_ITER_MEAN = zeros(NCases,2);
LDIPM_ITER_STD = zeros(NCases,2);
INFEAS_NUM = zeros(NCases,1);
for i = 1 : NCases
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(theFiles(i).folder, baseFileName);
    % fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    load(fullFileName); % Loads saveData structure
    if i == 1
        DATA = repmat(saveData,NCases,1); % initializes data structure
    else
        DATA(i) = saveData;
    end
    CG_ITER_MEAN(i,:) = mean(saveData.CG_ITERS)/1e3;
    CG_ITER_STD(i,:) = std(saveData.CG_ITERS)/1e3;
    LDIPM_ITER_MEAN(i,:) = mean(saveData.LDIPM_ITERS);
    LDIPM_ITER_STD(i,:) = std(saveData.LDIPM_ITERS);
    INFEAS_NUM(i) = sum(~saveData.FEASFLAG);
    SIZES(i,:) = [saveData.n saveData.m saveData.n_qp saveData.m_qp];
end

% Sort
[~,indVec] = sort(SIZES(:,3));
% indVec = indVec([2, 4, 6:end]);
CG_ITER_MEAN = CG_ITER_MEAN(indVec,:);
CG_ITER_STD = CG_ITER_STD(indVec,:);
LDIPM_ITER_MEAN = 2*LDIPM_ITER_MEAN(indVec,:);
LDIPM_ITER_STD = LDIPM_ITER_STD(indVec,:);
INFEAS_NUM = INFEAS_NUM(indVec,:);
SIZES = SIZES(indVec,:);
PERECENT = (CG_ITER_MEAN(:,1) - CG_ITER_MEAN(:,2))./CG_ITER_MEAN(:,1)*100
min(PERECENT)
max(PERECENT)
AVG = CG_ITER_MEAN*1e3./(LDIPM_ITER_MEAN)
