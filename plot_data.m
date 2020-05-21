clear all;

% read integer variables
int_vars = fopen('int_vars.txt','r');
formatSpec = '%d';
ints = fscanf(int_vars, formatSpec, [3 1]);
fclose(int_vars);
Nx = ints(1); Ny = ints(2); k = ints(3);

% read float variables
f_vars = fopen('f_vars.txt','r');
formatSpec = '%f';
floats = fscanf(f_vars, formatSpec, [7 1]);
fclose(f_vars);
dx = floats(1); dy = floats(2); SStime = floats(3);
Lx = floats(4); Ly = floats(5); T_max = floats(6);
T_min = floats(7);

x=zeros(1,Nx+2);y=zeros(1,Ny+2);            %Generate the plate
for i = 1:Nx+2
x(i) =(i-1)*dx;
end
for i = 1:Ny+2
y(i) =(i-1)*dy;
end

%% Plotting section
% %%%            -------------- Constant plot ----------------
% read Tss matrix
filename = 'tss.txt';
delimiterIn = '\t';
headerlinesIn = 0;
Tss = importdata(filename,delimiterIn,headerlinesIn);

figure(1);

hold on
title(sprintf('Temperature at steady state time : %i seconds ',round(SStime)))
surf(x,y,Tss)
cb=colorbar;
caxis([T_min T_max]);
view(90,-90);
xlim([0 Lx+dx]); xlabel('Length(m)');
ylim([0 Ly+dy]); ylabel('Width(m)');
zlim([T_min T_max]); zlabel('Temprature(C)');
drawnow
hold off

% %%%            -------------- Animated plot ----------------
% read T matrix
filename = 't.txt';
delimiterIn = '\t';
headerlinesIn = 0;
T = importdata(filename,delimiterIn,headerlinesIn);

figure (2);
j = 0; % iteration variable for movieVector

for i=1:20:k
    clf;
    hold on;
    title(sprintf('Temperature at time : %i seconds ',int32(0.6*i)));
    surf(x,y,T((Ny+2)*(i-1)+1:(Ny+2)*i,:));
    cb=colorbar;
    caxis([T_min T_max]);
    view(90,-90);
    xlim([0 Lx+dx]); xlabel('Length(m)');
    ylim([0 Ly+dy]); ylabel('Width(m)');
    zlim([T_min T_max]); zlabel('Temprature(C)');
    j = j+1;
    movieVector(j) = getframe(gcf);
end

all_valid = true;
flen = length(movieVector);
for K = 1 : flen
  if isempty(movieVector(K).cdata)
    all_valid = false;
    fprintf('Empty frame occurred at frame #%d of %d\n', K, flen);
  end
end
if ~all_valid
   error('Did not write movie because of empty frames')
end

myWriter = VideoWriter('plate');
myWriter.FrameRate = 20;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);

