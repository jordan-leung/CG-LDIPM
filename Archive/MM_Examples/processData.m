clc
clear all
close all

% Specify the folder where the files live.
myFolder = '/Users/jordan/Documents/GitHub/CG-LDIPM/MM_Examples/data';
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
N = length(theFiles);
FileNames = strings(N,1);
ProblemSize = zeros(N,1);
ConstraintSize = zeros(N,1);
ConditionNumber = zeros(N,1);
count = 0;
for k = 1 : N
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    % fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    load(fullFileName); % Loads variables A, p, l, m, n, q, r, u

    % Construct table elements
    if n <= 5000 && m <= 50000 && n >10 
        try chol(P);
            count = count + 1;
            FileNames(count) = convertCharsToStrings(baseFileName);
            ProblemSize(count) = n;
            ConstraintSize(count) = m;
        catch ME
            disp('Matrix is not symmetric positive definite')
        end
    end
    k
end

% Delete extra rows
FileNames = FileNames(1:count);
ProblemSize = ProblemSize(1:count);
ConstraintSize = ConstraintSize(1:count);
ConditionNumber = ConditionNumber(1:count);

TABLE = table(FileNames,ProblemSize,ConstraintSize,ConditionNumber);
save('TABLE','TABLE')