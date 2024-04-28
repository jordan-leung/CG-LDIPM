clc
clear all
close all

% Cases, stored as the size (n,m,targetConditionNumber)
n = 100;
m = 500;
condNum = 100;

% Number of samples
NSamples = 50;

% Zeros
A = zeros(m,n);

for j = 1:NSamples
    if j == 1 % initialize structure
        % Randomly sample A
        for ii = 1:m
            p_i = normrnd(zeros(1,n),1);
            a_i = p_i/norm(p_i);
            A(ii,:) = a_i;
        end

        % Get hh as a random NxN matrix
        r = n;
        R = zeros(r,n);
        for ii = 1:r
            p_i = normrnd(zeros(1,n),1);
            r_i = p_i/norm(p_i);
            R(ii,:) = r_i;
        end
        HPrime = R'*R;
        [u, s, v] = svd(HPrime);
        s = diag(s);           % s is vector
        s = s(1)*( 1-((condNum-1)/condNum)*(s(1)-s)/(s(1)-s(end))) ;
        s = diag(s);           % back to matrix
        H = u*s*v';
        H = 1/2*(H' + H);

        % Get the vectors
        x = normrnd(zeros(n,1),1);
        w1 = normrnd(zeros(m,1),1);
        w2 = normrnd(zeros(m,1),1);
        s = ones(m,1) + 0.1*abs(w1);
        lambda = ones(m,1) + 0.1*abs(w2);
        b = s - A*x;
        c = A'*lambda - H*x;

        % Form the structure for initialization
        dataInit.W = H;
        dataInit.c = c;
        dataInit.A = A;
        dataInit.b = b;

        % Now, form the NCases x NSubcases structure array
        dataStruct = repmat(dataInit,NSamples,1); % the (1,1) case is auto stored
    else
        % Randomly sample A
        A = zeros(m,n);
        for ii = 1:m
            p_i = normrnd(zeros(1,n),1);
            a_i = p_i/norm(p_i);
            A(ii,:) = a_i;
        end

        % Get hh as a random NxN matrix
        r = n;
        R = zeros(r,n);
        for ii = 1:r
            p_i = normrnd(zeros(1,n),1);
            r_i = p_i/norm(p_i);
            R(ii,:) = r_i;
        end
        HPrime = R'*R;
        [u, s, v] = svd(HPrime);
        s = diag(s);           % s is vector
        s = s(1)*( 1-((condNum-1)/condNum)*(s(1)-s)/(s(1)-s(end))) ;
        s = diag(s);           % back to matrix
        H = u*s*v';
        H = 1/2*(H' + H);

        % Get the vectors
        x = normrnd(zeros(n,1),1);
        w1 = normrnd(zeros(m,1),1);
        w2 = normrnd(zeros(m,1),1);
        s = ones(m,1) + 0.1*abs(w1);
        lambda = ones(m,1) + 0.1*abs(w2);
        b = s - A*x;
        c = A'*lambda - H*x;

        % Form the structure for initialization
        dataStruct(j).W = H;
        dataStruct(j).c = c;
        dataStruct(j).A = A;
        dataStruct(j).b = b;
    end
end
