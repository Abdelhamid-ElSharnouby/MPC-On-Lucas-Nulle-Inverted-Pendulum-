function [] = observability(A,cT)
%OBSERVABILITY checks whether the system is observable 
%   observability uses the function OBSV to compute the 
%   observability matrix. The system is observable 
%   if the observability matrix has full rank n. 
n = size(A,1);
% observability
if isequal(rank(obsv(A,cT)),n);
    disp('Your system is completely observable!')
else
    disp('Your system is not observable!')
end


end
