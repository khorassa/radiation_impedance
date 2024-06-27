% Numerical Integration script by Eric Dew for the radiation impedance of
% rectangular CMUTs
%   Q is an array of all the aspect ratios you want to include
%   kax is the range of ka values (frequency values) to integrate over (for the 2D case)
%   kax1D is the same for the 1D case which is slower and may benefit from fewer points.
%   Simulate1D and 2D are options to toggle off simulating one of the cases.
% Note that the output of this script is dimensionless, in other words, it
% computes Z/(rho * c * 4ab); where Z is the radiation impedance, rho is
% the medium density, c is the speed of sound, and 4ab is the total area of the membrane
% The results are saved into a directory named Rad_Imp
mkdir('Rad_Imp\OneDZ');
mkdir('Rad_Imp\TwoDZ');

Q = [1 4 10 25]; % aspect ratios to simulate

kaxInc = 0.025; % increment in kax in Log scale
kax = 10.^(-2:kaxInc:log10(20));
kaxInc1D = 0.1; % increment in kax1D in Log scale
kax1D = 10.^(-2:kaxInc1D:1); % going up to 10 since 1D is slower
xMin = 0.1 % the threshold under which the Taylor expansion is used instead
% of the original expression, because the original expression suffers from precision errors

% Set to 1 to compute and 0 to not
Simulate2D = 1;
Simulate1D = 1;

if Simulate2D
    % define a function handle SP that calls Spw - custom piecewise function with the infinities dealt with and a taylor expn at 0
    SP = @Spw;
    % define function handle Sp that calls SP with argument x - for convenience
    Sp =@(x)SP(x,xMin);
    fprintf('\nSimulating 2D Radiation Impedance Results')
    for jj = 1:length(Q)
        q = Q(jj); % aspect ratio: b/a
        fprintf('\nSimulating Aspect ratio: q = %.3f. Number %i out of %i.',q,jj,length(Q))
        RR = zeros(length(kax),1);
        XX = zeros(length(kax),1);
        for ii = 1:length(kax)
            kb = q*kax(ii);
            INTFUN = @(t,phi)   Sp(kax(ii) * t .* cos(phi)) .^2  .* Sp(kb * t .* sin(phi)) .^2;
            INTFUNR = @(t,phi) INTFUN(t,phi) .* t ./ sqrt(1 - t.^2); %for R / 0<=t<=1
            INTFUNX = @(t,phi) INTFUN(t,phi) .* t ./ sqrt(t.^2 - 1); %for X / 1<=t<inf
            % due to symmetry, we can integrate from 0 to pi/2 and multiply by 4
            Rtest = 4*integral2(INTFUNR,0,1,0,pi/2);
            Xtest = 4*integral2(INTFUNX,1,inf,0,pi/2);
            RR(ii) = kax(ii)*kb * 1/(16 * pi^2) * (315/128)^2 *  Rtest;
            XX(ii) = -1* kax(ii)*kb * 1/(16 * pi^2) * (315/128)^2 *  Xtest;
            fprintf('\nFinished Integration of ka = %.3f. ',kax(ii))
            fprintf('%i out of %i completed',ii,length(kax))
        end
        writematrix(RR,strcat('Rad_Imp\TwoDZ\',sprintf('RRqis%d.txt',q)));
        writematrix(kax,strcat('Rad_Imp\TwoDZ\',sprintf('KAXqis%d.txt',q)));
        writematrix(XX,strcat('Rad_Imp\TwoDZ\',sprintf('XXqis%d.txt',q)));
        fprintf('\n')
    end
end

if Simulate1D
    % define a function handle SP that calls Spw - custom piecewise function with the infinities dealt with and a taylor expn at 0
    SP = @Spw;
    % define function handle Sp that calls SP with argument x - for convenience
    Sp =@(x)SP(x,xMin);
    % 1D Version:
    CS = @cSinc; %handle for sinc:
    
    fprintf('\nSimulating 1D Radiation Impedance Results')
    parfor jj = 1:length(Q)
        q = Q(jj); %aspect ratio: b/a
        fprintf('\nSimulating Aspect ratio: q = %.3f. Number %i out of %i.',q,jj,length(Q))
        RR1D = zeros(length(kax1D),1);
        XX1D = zeros(length(kax1D),1);
        for ii = 1:length(kax1D)
            kb = q*kax1D(ii);
            INTFUN1D = @(t,phi)   Sp(kax1D(ii) * t .* cos(phi)) .^2  .* CS(kb * t .* sin(phi)) .^2;
            INTFUNR1D = @(t,phi) INTFUN1D(t,phi) .* t ./ sqrt(1 - t.^2); %for R / 0<=t<=1
            INTFUNX1D = @(t,phi) INTFUN1D(t,phi) .* t ./ sqrt(t.^2 - 1); %for X / 1<=t<inf
            % due to symmetry, we can integrate from 0 to pi/2 and multiply by 4
            Rtest = 4*integral2(INTFUNR1D,0,1,0,pi/2);
            Xtest = 4*integral2(INTFUNX1D,1,inf,0,pi/2);
            RR1D(ii) = kax1D(ii)*kb * 1/(4 * pi^2) * (315/128) *  Rtest;
            XX1D(ii) = -1* kax1D(ii)*kb * 1/(4 * pi^2) * (315/128) *  Xtest;
            fprintf('\nFinished Integration of ka = %.3f. (q = %i)',kax1D(ii),q)
            fprintf('%i out of %i completed. (q = %i)',ii,length(kax1D),q) 
        end
        writematrix(RR1D,strcat('Rad_Imp\OneDZ\',sprintf('RRqis%d.txt',q)));
        writematrix(kax1D,strcat('Rad_Imp\OneDZ\',sprintf('KAXqis%d.txt',q)));
        writematrix(XX1D,strcat('Rad_Imp\OneDZ\',sprintf('XXqis%d.txt',q)));

    end

end






function s =  Spw (x, xMin)
    %piecewise version of S(x)
    S =@(x) -16 * (3*x.* cos(x) + (x.^2 -3).*sin(x)  ) ./ x.^5;
    S_tay = @(x) (2.*x.^4) /945 - (8*x.^2)/105 + 16/15;
    %global max of S is 16/15, so anything that appears weird from numerical errors has a max value of that.
    SF = S(x);
    SF(isnan(SF)) = 16/15; %just set this to something finite since the taylor series should be used
    %Could implement if isnan(SF) and ~isnan(ST) set it to zero
    SF(isinf(x)) = 0; %the limit at inf is 0
    
    ST = S_tay(x);
    ST(isinf(x)) = 0; %enforce limit on taylor series else we also get this problem of adding 0*NaN
    
    %if x is NaN: eg: S( t-> inf * sin(0)). - realistically this will go to zero, but technically S(0*inf) is undefined.
    %max of the function is 16/15 so we can set it to that, in that case it should always be multiplied by zero in the overall integral anyways, so that edge case
    %should never actually give a non-zero value. But I'll set it to the global max of the function to be as conservative as possible about not zeroing things unneccsarily.
    ST(isnan(x)) = 0; %zero it in the taylor expansion because for some reason I am still summing the 2 pieces instead of using an if statement.
    SF(isnan(x)) = 16/15; %upper bound on the value but the overall function should be zero in this case
    
    %make sure only the correct expression is used
    ST(abs(x)>=xMin) = 0;
    SF(abs(x)<xMin) = 0;
    %silly way to do piecewise by just adding them and doing 0* the other one.
    s = (abs(x) < xMin ) .* ST + (abs(x) >= xMin) .* SF;
end


function s =  cSinc (x)
    %custom sinc function;
    CSinc = @(x) sinc(x/pi); %get rid of MATLAB's sinc(x) = sin(x*pi)/(x*pi) convention.
    SF = CSinc(x);
    %global max of sinc is 1, so anything that appears weird from numerical errors has a max value of that.
    SF(isinf(x)) = 0; %the limit at inf is 0
    %if x is NaN: eg: CS( t-> inf * sin(0)). - realistically this will go to zero, but technically CS(0*inf) is undefined.
    %max of the function is 1 so we can set it to that, in that case it should always be multiplied by zero in the overall integral anyways, so that edge case
    %should never actually give a non-zero value. But I'll set it to the global max of the function to be as conservative as possible about not zeroing things unneccsarily.
    SF(isnan(SF)) = 1; %upper bound on the value but the overall function should be zero in this case
    
    s = SF;
end
