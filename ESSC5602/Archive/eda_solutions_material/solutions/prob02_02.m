% Problem 2.1

% read data from file
D=load('brf_temp.txt');
t=D(:,1);
d=D(:,2);
N = length(t);
Dt = t(2)-t(1);

% initialize largest change
largest=0.0;
timeoflargest=t(1);

% The loops computes the derivative dddt at each time, i.
% Only good data are used, the others are skipped over.
for i=[1:N-1]
    
    % check that no times are missing between
    % t(i) and t(i+1).  The time difference between
    % samples i and i+1 should be about Dt.
    thisDt = (t(i+1)-t(i));
    if( abs((thisDt-Dt)/Dt) > 0.1  )
        continue;
    end
    
    % check if d(i) is a bad sample
    if( (d(i)<-40) || (d(i)>38) || (d(i)==0) )
        continue;
    end
    
    % check if d(i+1) is a bad sample
    if( (d(i+1)<-40) || (d(i+1)>38) || (d(i+1)==0) )
        continue;
    end 
    
    % compute derivative, units of K/day
    dddt = (d(i+1)-d(i)) / thisDt;
    
    % update largest if necessary
    if( abs(dddt) >= largest )
        largest = abs(dddt);
        timeoflargest=t(i);
    end
    
end

disp(sprintf('largest hourly change is %f K/hr occurs at time %f days',largest/24,timeoflargest));




