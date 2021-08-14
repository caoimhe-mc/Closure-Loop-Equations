function [x,y,N,err] = closureThreeBar(a,b,c)
 %Calculates two unknown angles based on a loop closure equation 
 %for a three-bar linkage system where all the lengths of the linkage system 
 %are known
 %
 % Inputs:
 % a,b,c (real scalar, positive)lengths of linkages
 %units for the lengths can vary, however each length must be 
 %expressed in the same units
 % Outputs:
 % x the angle between links a and c in radians
 % y the angle between links a and b in radians
 % N is the number of iterations completed before terminating.
 % err is the error  
 
 %Version 1: created 15/03/2021. Author: C.McCann
 
 %determines whether the input is a positive,real and a scalar number and an error 
 %is diplayed if the input does not meet the criteria
for m = [a b c]
    if (~isscalar(m)) || (~isreal(m)) || m <= 0
        error('Input is invalid')
    end
end
 
x = pi/4; %inital guess for angle between link a and c
y = pi/4; %inital guess for angle between link a and b
k = [x,y];
tolerance = 10e-7; %the limit to determine when to terminate
iteration = 20;    %defines the number of iterations if the tolerance 
                   %is never reached                   

 for N = 1: iteration+1

     %terminates the loop if the values do not converge within the limit
     if N == iteration + 1  
         error('Iteration limit reached. Iteration did not converge.')
         break
     end

     %updates the two scalar equations determined from the vector eqaution
     f = c*cos(k(1)) + b*cos(k(2)) -a; 
     g = c*sin(k(1)) - b*sin(k(2));        
     func = [f;g]; %a column vector of the two functions

     %terminates the loop if the absolute values are within the tolerance value    
     if abs(norm(func)) < tolerance
          break
     end

     dfx = -c*sin(k(1)); %the derivative of the first function with respect to x
     dfy = -b*sin(k(2)); %the derivative of the first function with respect to y
     dgx = c*cos(k(1));  %the derivative of the second function with respect to x
     dgy = -b*cos(k(2)); %the derivative of the second function with respect to y
     Jacobian = [dfx dfy; dgx dgy]; %the Jacobian is a 2x2 matirx
     k = k - Jacobian\func; %updates the values for x and y

 end

 x = k(1);
 y = k(2);
 err = abs(norm(f));

end