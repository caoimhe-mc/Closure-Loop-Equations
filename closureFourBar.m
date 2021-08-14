function [x,y,N,err] = closureFourBar(a,b,c,d,angle)
 %Calculates two unknown angles in a four-bar linkage based on a
 %closure loop equation where all the lengths of the linkage system are 
 %known
 %
 % Inputs:
 %a,b,c,d (real scalar, positive)lengths of linkages 
 %units for the lengths can vary, however each length must be 
 %expressed in the same units
 %angle is the known angle required in radians
 % Outputs:
 % x the angle between links d and a in radians
 % y the angle between links c and the horizontal in radians
 % N is the number of iterations completed before terminating.
 % err is the error  
 
 %Version 1: created 15/03/2021. Author: C.McCann
 
 %determines whether the input is a positive,real and a scalar number and an error 
 %is diplayed if the input does not meet the criteria
for m = [a b c d]
    if (~isscalar(m)) || (~isreal(m)) || m <= 0
        error('Input is invalid')
    end
end

x = 5*pi/9;
y = 19*pi/18;
k = [x,y];
tolerance = 10e-8; %the limit to determine when to terminate
iteration = 20;    %defines the number of iterations if the tolerance 
                   %is never reached                   

 for N = 1: iteration+1
     
     %terminates the loop if the values do not converge within the limit
     if N == iteration + 1  
         error('Iteration limit reached. Iteration did not converge.')
         break
     end
    
     %updates the two scalar equations determined from the vector eqaution
     f = a*cos(k(1)) - b*cos(k(2)) - c*cos(angle)-d;
     g = a*sin(k(1)) - b*sin(k(2)) - c*sin(angle);
     func = [f;g]; %a column vector of the two functions
     
     if abs(norm(func)) < tolerance
          break
     end
     
     dfx = -a*sin(k(1)); %the derivative of the first function with respect to x
     dfy = b*sin(k(2));  %the derivative of the first function with respect to y
     dgx = a*cos(k(1));  %the derivative of the second function with respect to x
     dgy = -b*cos(k(2)); %the derivative of the second function with respect to y
     Jacobian = [dfx dfy; dgx dgy]; %the Jacobian is a 2x2 matirx 
     k = k - Jacobian\func; %updates the values for x and y
     
 end
 
 x = k(1);
 y = k(2);
 err = abs(norm(f));
 
end