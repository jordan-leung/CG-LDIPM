clc
clear all
close all

n = 50;
m = 150;
condTarget = 100;

% Randomly sample A
A = zeros(m,n);
for i = 1:m
    p_i = normrnd(zeros(1,n),1);
    a_i = p_i/norm(p_i);
    A(i,:) = a_i;
end

% Get hh as a random NxN matrix
r = n;
R = zeros(r,n);
for i = 1:r
    p_i = normrnd(zeros(1,n),1);
    r_i = p_i/norm(p_i);
    R(i,:) = r_i;
end
HPrime = R'*R;
[u, s, v] = svd(HPrime);
s = diag(s);           % s is vector
s = s(1)*( 1-((condTarget-1)/condTarget)*(s(1)-s)/(s(1)-s(end))) ;
s = diag(s);           % back to matrix
H = u*s*v';
H = 1/2*(H' + H);


x = normrnd(zeros(n,1),1);
w1 = normrnd(zeros(m,1),1);
w2 = normrnd(zeros(m,1),1);
s = ones(m,1) + 0.1*abs(w1);
lambda = ones(m,1) + 0.1*abs(2);
b = s - A*x;
c = A'*lambda - H*x;

save('randData2','H','c','A','b')