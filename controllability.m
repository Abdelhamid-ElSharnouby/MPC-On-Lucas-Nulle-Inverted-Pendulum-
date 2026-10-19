function [] = controllability(A,b)
%CONTROLLABILITY checks whether the system is controllable 
%   controbility uses the function CTRB to compute the 
%   controllability matrix. The system is controllable 
%   if the controllability matrix has full rank n. 
n = size(A,1);
% controllability
if isequal(rank(ctrb(A,b)),n);
    disp('Your system is completely controllable!')
else
    disp('Your system is not controllable!')
end
end
