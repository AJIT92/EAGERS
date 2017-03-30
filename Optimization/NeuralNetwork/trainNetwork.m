function [Net,sqrerror] = trainNetwork(Net,desiredOut, inputs)
%this does forward propagation for a one layer network for a set of
%generators and a demand
%inputs: network, desiredOut: desired network output when using forward
%funnction in the form of a vertical vector, inputs: matrix of inputs of
%size inputlength x number of outputs
%inputs in order: ub, lb, f, H for each generator, demand, $/kWgrid

[sqrerror,dedW,dedb] = finderror(Net,inputs,desiredOut);%find the error and the gradient of the error

%% perturb weights to check that derrordW is found correctly before trying to train using derrordW
numgrad = zeros(size(Net.Wlayer1));
perturb = numgrad;
Winitial = Net.Wlayer1;
netup = Net;
netdwn = Net;
for j = 1:1:length(numgrad(1,:))
    for i = 1:1:length(numgrad(:,1))
        perturb(i,j) = .01*1/Net.nodeconst;
        netup.Wlayer1 = Winitial+perturb;
        netdwn.Wlayer1 = Winitial-perturb;
        [sqrerup,~,~] = finderror(netup, inputs, desiredOut);
        [sqrerdwn,~,~] = finderror(netdwn, inputs, desiredOut);
        numgrad(i,j) = (sum(sum(sqrerup-sqrerdwn)))/(2*perturb(i,j));
        perturb(i,j) = 0;
    end
end
gradientdiff = norm(dedW-numgrad)/norm(dedW+numgrad);
if gradientdiff>1e-6
    disp('Warning, gradient function may be incorrect. Training may be inaccurate.')
end

%% now that you have already determined that derrordW is calculated correctly
%% use BFGS technique to train using derrordW
% if Net.classify
    tolerance = .0001.*ones(1,length(sqrerror(1,:)));
% else
%     tolerance = 1.*ones(1,length(sqrerror(1,:)));
% end
ddeddW = -dedW;%find second derivative of error, 
% initialize an approximation of the Hessian matrix = d^2f/(dx_i dx_j)
warning('off','all')%prevent print of warning as Hessian gets close to singular
for i = 1:1:length(dedW(1,:)) %each set for each node output must be trained individually
%     H = eye(length(dedW(:,1)));% H*trainingstep = - dedW
    sqrerror1 = sqrerror(i,:);
    dedW1 = dedW(:,i);
    dedb1 = dedb(:,i);
    iterations = 0;
    laststep = zeros(size(dedW1));
    lastbstep = zeros(size(dedb1));
    a = 100;
    failed = false;
%     a = 100/(100^Net.classify);%start at 100 for numeric, start at 1 for classification
    while nnz(sqrerror1>tolerance)>0 %keep training until you get the desired output
        iterations = iterations+1;
        trainingstep = dedW1;
        bstep = dedb1;%training step for bias
        %trainingstep = H\dedW1;%inv(H)*dedW1; %obtain training step direction
        %perform line search to find acceptable scalar for step size
        %must find scalar that minimizes error from Wlayer1+scalar*trainingstep
        %solve 0 = sqrerror(Wlayer1+scalar*trainingstep), or require sufficient
        %decrease in error
        if nnz(isnan(trainingstep))>0 || nnz(isinf(trainingstep))>0 || nnz(isnan(bstep))>0 || nnz(isinf(bstep))>0%nnz(isnan(dedWnew))>0 || nnz(isnan(sqrerrornew))>0 || 
%             disp('warning, training returns NaN')
%             trainingstep(trainingstep==-inf) = dedW1(trainingstep==-inf)*(-10e10);
%             trainingstep(isinf(trainingstep)) = dedW1(isinf(trainingstep))*10e10;
            trainingstep(trainingstep==-inf) = -10e10;
            trainingstep(isinf(trainingstep)) = 10e10;
            trainingstep(isnan(trainingstep)) = 0;
            bstep(bstep==-inf) = -10e10;
            bstep(isinf(bstep)) = 10e10;
            bstep(isnan(bstep)) = 0;
        end
        scalar = a;%/max(max(abs(trainingstep)));
        if Net.classify
            momentum = .1;%.3 is too high, .1 is too low, .2 does well for test2E_1BS
        else momentum = 0.5;
        end
        step = -trainingstep.*scalar+laststep.*momentum;
        bstep = -bstep.*scalar/100+lastbstep.*momentum/100;
        Net.Wlayer1(:,i) = Net.Wlayer1(:,i)+step;
        Net.blayer1(:,i) = Net.blayer1(:,i)+bstep;
        [sqrerrornew, dedWnew, dedbnew] = finderror(Net, inputs, desiredOut);
        sqrerrornew1 = sqrerrornew(i,:);
        dedWnew1 = dedWnew(:,i);
        dedbnew1 = dedbnew(:,i);
        if sum(sum(abs(sqrerrornew1)))>=sum(sum(abs(sqrerror1))) || nnz(abs(dedWnew1)>1e-30)==0 || nnz(isinf(sqrerrornew1))>0%|| nnz(isnan(dedWnew))>0%if the error gets worse or you have reached a flat point
            if abs(a) <1e-12 %try the other direction
                a = -a*1e15/(100^Net.classify);
            else
                a = a/10;
            end
            Net.Wlayer1(:,i) = Net.Wlayer1(:,i)-step;
            Net.blayer1(:,i) = Net.blayer1(:,i)-bstep;
            laststep = zeros(size(laststep));
            lastbstep = zeros(size(lastbstep));
        else
%             yk = dedWnew1 - dedW1;%this is the change in dedW from the step
            laststep = step;
            lastbstep = bstep;
%             if nnz(abs(dedWnew1)>1e-30)==0
%                 disp('reached flat gradient, exiting loop');
%                 sqrerror(i,:) = sqrerrornew1;
%                 break
%             else
                if nnz(sqrerrornew1>tolerance)==0
                    sqrerror(i,:) = sqrerrornew1;
                    disp('below tolerance');
                end
                sqrerror1 = sqrerrornew1;
                dedW1 = dedWnew1;
                dedb1 = dedbnew1;
                %find the new Hessian
%                 if nnz(yk)==0 %reached linear section, move really fast
%                     a = a*1e5;
%                 end
%                 H = H + yk*yk'/(yk'*step) - H*(step*step')*H/(step'*H*step);%H_k+1 = H_k + U_k + V_k
%                 if nnz(isinf(H))>0 || nnz(isnan(H))>0 %|| rcond(H)<1e-15%(det(H)<1e-16 && det(H)>-1e-16)
% %                     disp('warning: Hessian is NaN or close to singular. return to identity matrix.')
%                     H = eye(length(step));                    
%                 end
%             end
        end
        if iterations>1e+4*length(desiredOut(1,:))/100%10^4 iterations per 100 timesteps
            disp('not converging after 10^4 iterations, exiting loop');
            if Net.classify
                if failed ==false && sum(sum(sqrerror))>length(sqrerror)*0.025%wrong for more than 5 percent of timesteps
                    Net.Wlayer1(:,i) = rand(length(Net.Wlayer1(:,i)),1);%start with new initial values
                    iterations = 0;
                    failed = true;
                else
                    sqrerror(i,:) = sqrerrornew1;
                    break
                end

            else
                sqrerror(i,:) = sqrerrornew1;
                break
            end
        end
    end
end



function [cost,derrordW,derrordb] = finderror(Net,inputs,desiredOut)
NetOut = forward(Net,inputs);
error = (desiredOut-NetOut);%all errors
cost = error.^2.*0.5;%/length(inputs(:,1)) + (Net.lambda/2)*sum(Net.Wlayer1)'.^2*ones(1,length(desiredOut(1,:)));% normalized(1/2error^2) + 1/2penalty for complex model
%keep model simple using lambda to prevent over fitting
if Net.classify %if it has a sigmoid function
    %cost = -1/length(error(1,:))*(desiredOut.*log(NetOut)+(1-desiredOut).*log(1-NetOut));
    %derrordW = (-error.*Net.nodeconst.*NetOut./(1+exp((inputs*Net.Wlayer1)'.*Net.nodeconst))*inputs)';
    %use cross error to prevent learning slowdown with sigmoid functions
    derrordW = -1/length(error(1,:))*(error*inputs)';
else%if no activation function
    derrordW = (-error*inputs)';% + Net.lambda*Net.Wlayer1;%no activation function so this is just the error
end
derrordb = -1/length(error(1,:))*sum(error,2)';
