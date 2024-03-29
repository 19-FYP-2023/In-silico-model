function absSum = model(parameters)
%MODEL Summary of this function goes here
%   Detailed explanation goes here
%   parameters
%   1 - diameter

    nx = 200;
    ny = 200;
    nz = 200;
    x = linspace(-0.3, 0.3, nx);
    y = linspace(-0.3, 0.3, ny);
    z = linspace(-0.8, 0, nz);
    % x = 1:1:nx;
    % y = 1:1:ny;
    % z = 1:1:nz;
    [X,Y,Z] = ndgrid(single(x),single(y),single(-z));

    [mx, my, mz] = meshgrid(1:1:100);
    [M,A] = geometryDefinition(X,Y,Z,parameters);

    d = 0.3;
    h = 0.4;
    a = (4*h)/(d^2);
    e = a.*(x - (d/2)).*(x + (d./2));
    e = min(0,e);
    
    ro = 0.01;
    norm_amp = 0.02;

    px = (x>-d/2 & x<d/2);
    absSum = 0;

%     absarr = [];


    for i= 1:length(x)
        if px(i) == 1
           mu = [0 e(i)];
           Sigma = [ro.*abs(e(i)) 0;0 ro.*abs(e(i))];
           [YM,ZM] = meshgrid(y,z);
           YZ = [YM(:) ZM(:)];
           norm = norm_amp.*mvnpdf(YZ,mu,Sigma);
           norm = reshape(norm,length(z),length(y))';

%            TM = M(:,i,:);
%            TM = reshape(TM,length(z),length(y));

           TA = A(:,i,:);
           TA = reshape(TA,length(z),length(y));

           thresh = 0.05;

           temp = TA.*norm.*(norm>thresh);
           absSum = absSum + sum(temp(:));
%            absarr(end+1,1) =  sum(TA.*norm.*(norm>thresh),"all");
%            M(:,i,:) = norm.*(norm>thresh) + TM.*(norm<=thresh);
%            max_val = max(max(norm));

        end
    end


end

