function [ ] = robustness( A )

%ROBUSTNESS checks the robustness of A state-space model
%   The function calculate the eigenvalues of the state-space model.
%   This will determine whether it is a stable system or not. 
pole = eig(A);
imagPol = find(imag(pole)); 
center = find(real(pole)==0);

    if max(real(pole)) > 0
        disp('Your system is unstable!')
    elseif max(real(pole)) < 0
        disp('Your system is stable!')
    elseif length(center)>2
        disp('Your system is unstable!')
    elseif length(center)-length(intersect(center,imagPol))>1  
        disp('Your system is unstable!')
    elseif length(intersect(center,imagPol))==2
        disp('Your system is critically stable leads harmonic continuous oscillations!')
    else
        disp('Your system is critically stable!')
    end
   
end


        
    