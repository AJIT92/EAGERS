function interpData = interpolateData(Tfreq,Xi,Xf,Variability)
global TestData
S = fields(TestData.Demand);
r = ceil((TestData.Timestamp(2)-TestData.Timestamp(1))*24*3600/Tfreq);% # of points that must be created between datum
lengthTD = length(TestData.Timestamp);
%% create noise vector = 10*# of sub steps
if r>0 
    z = randn(round(100*r),1);
    N = zeros(round(100*r),1);
    N(1) = 0;
    N(2) = (rand(1,1)-.5)*Variability; %the value .04 makes the final noise signal peaks = variability with this strategy.
    b= sign(N(2)-N(1));
    c = N(2);
    for n = 3:length(N)
        %if the noise is increasing the probability is such that the noise should continue to increase, but as the noise approaches the peak noise magnitude the probability of switching increases.
        if (c>0 && N(n-1)>0) || (c<0 && N(n-1)<0)
            a = 2*(Variability - abs(N(n-1)))/Variability; %constant 2 = 97.7% probability load changes in same direction, 1.5 = 93.3%, 1 = 84.4%, .5 = 69.1%
        else a = 2;
        end
        if abs(z(n))>.68 %only changes value 50% of the time
             c = b*(z(n)+a); % c is positive if abs(Noise)is increasing, negative if decreasing
             b = sign(c);
             N(n) = N(n-1)+c*Variability; 
        else N(n) = N(n-1);
        end
    end
    Noise = N;% scaled noise to a signal with average magnitude of 1
end

%% take data exactly, or interpolate if necessary
if r==1
    interpData.Timestamp=TestData.Timestamp(Xi:Xf);
    interpData.Temperature=TestData.Temperature(Xi:Xf);
    for i = 1:1:length(S)
        interpData.Demand.(char(S(i)))=TestData.Demand.(char(S(i)))(Xi:Xf);
    end
elseif r<1 %extra datum, use nearest timepoint
    NumSteps = (TestData.Timestamp(Xf)-TestData.Timestamp(Xi))*24*3600/Tfreq+1;
    interpData.Timestamp =[];
    interpData.Temperature = zeros(NumSteps,1);
    x1 = Xi;
    for i = 1:1:NumSteps
        x2=nnz(TestData.Timestamp<(TestData.Timestamp(Xi)+i*Tfreq/(24*3600)));
        interpData.Temperature(i) = mean(TestData.Temperature(x1:x2));
        interpData.Timestamp(i) = TestData.Timestamp(x1);
        for j = 1:1:length(S)
            interpData.Demand.(char(S(j)))(i) = mean(TestData.Demand.(char(S(j)))(x1:x2));
        end
        x1=x2+1;
    end
elseif r>1 %interpolate between timesteps
    NumSteps = floor((TestData.Timestamp(Xf)-TestData.Timestamp(Xi))*24*3600/Tfreq+1);
    interpData.Timestamp =zeros(NumSteps,1);
    interpData.Temperature = zeros(NumSteps,1);
    for j = 1:1:length(S)
        interpData.Demand.(char(S(j))) = zeros(NumSteps,1);
    end
    nn = length(Noise);
    for j = 1:1:length(S)
        x1 = Xi;
        i=0;
        xi = ceil((length(Noise)-1000)*rand); %random starting point in extra long noise signal
        prev = TestData.Demand.(char(S(j)))(x1);
        x2 = x1+1;
        next = TestData.Demand.(char(S(j)))(x2);
        while i<NumSteps
            i=i+1;
            ni = mod(xi+i-1,nn)+1;
            r0 = (mod(i-1,r))/r;
            if j==1
                interpData.Timestamp(i) = TestData.Timestamp(Xi)+(i-1)*Tfreq/(24*3600);
                interpData.Temperature(i) = TestData.Temperature(x1)*(1+Noise(ni));
            end
            if TestData.Timestamp(x2)<=interpData.Timestamp(i)
                x1 = x1+1;
                prev = next;
                x2 = min(x1+1,lengthTD);
                next = TestData.Demand.(char(S(j)))(x2);
            end
            interpData.Demand.(char(S(j)))(i) = ((1-r0)*prev+r0*next)*(1+Noise(ni));
        end
    end
end