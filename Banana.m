fs = 20;  % Sampling rate (in Hz)
t = 0:1/fs:1;  % Time vector (replace 10 with the duration of your PPG signal)
ppg_data = (sin(2*pi*1*(t+0.85)) + 0.5*sin(2*pi*2*(t+0.85)) + 0.2*sin(2*pi*3*(t+0.85)))/1.4

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
[M,A] = geometryDefinition(X,Y,Z,{ppg_data(1)*0.05});

d = 0.3;
a = 18;
e = a.*(x - (d/2)).*(x + (d./2));
e = min(0,e);


px = (x>-d/2 & x<d/2);
absSum = 0;

absarr = [];


for i= 1:length(x)
    if px(i) == 1
       mu = [0 e(i)];
       Sigma = [0.01.*abs(e(i)) 0;0 0.01.*abs(e(i))];
       [YM,ZM] = meshgrid(y,z);
       YZ = [YM(:) ZM(:)];
       norm = 0.02.*mvnpdf(YZ,mu,Sigma);
       norm = reshape(norm,length(z),length(y))';
      
       TM = M(:,i,:);
       TM = reshape(TM,length(z),length(y));

       TA = A(:,i,:);
       TA = reshape(TA,length(z),length(y));
       
       thresh = 0.05;
       
       
       absSum = absSum + sum(TA.*norm.*(norm>thresh),"all");
       absarr(end+1,1) =  sum(TA.*norm.*(norm>thresh),"all");
       M(:,i,:) = norm.*(norm>thresh) + TM.*(norm<=thresh);
       
%        disp("ei");
%        disp(i);
%        disp("mean");
       max_val = max(max(norm));
%        [maxx,maxy] = find(norm==max_val);
%        disp(max_val);
    end
    
end

disp("absSum");
disp(absSum);


figure
axis 'equal';
h = slice(x,y,z,M,[],0,[]);
set(h,'edgeColor','none');
hold 'on';

hz = slice(x,y,z,M,0,[],[]);
set(hz,'edgeColor','none');
hold 'on';
% 
% hy = slice(x,y,z,M,[],[],-0.4);
% set(hy,'edgeColor','none');
% hold 'on';

colormap(turbo)

plot3(x,zeros(ny),e);
hold 'on';

% for i= 1:length(x)
%     y0 = 0;
%     z0 = e(i);
%     r = abs(min(0,5.*(x(i) - (d/2)).*(x(i) +(d./2))));
%     theta = linspace(0,2*pi,100);
%     %plot3(x(i).*ones(100),y0 + r*cos(theta),z0 + r*sin(theta),'-')
%     hold 'on';
%     
% end


% px = linspace(-d/2, d/2, 102);
% 
% for i= 1:length(px)
%    vx = interp1(x,e,px(i));
%    scatter3(px(i),0,vx,'filled')
%    hold 'on';
% end


b = a+5;
e = b.*(x - (d/2)).*(x + (d./2));
e = min(0,e);
plot3(x,zeros(ny),e);
hold 'on';

b = a-5;
e = b.*(x - (d/2)).*(x + (d./2));
e = min(0,e);
plot3(x,zeros(ny),e);
hold 'on';

figure 
plot(absarr);

