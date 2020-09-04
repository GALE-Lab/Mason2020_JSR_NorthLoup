% Calculation of deformation and other parameters in bedforms mapped in the
% North Loup river NE.

% The following script accompanies the publication entitled 
% Subaqueous dune field pattern evolution and interactions: North Loup 
% River, Nebraska, USA" by Mason et al. currently in review with the
% Journal of Sedimentary Research. For questions please contact
% corresponding author Mackenzie Day at daym@epss.ucla.edu 

% This script takes as an input an XLS file with columns defining dune
% label, and x and y points on that dune, and sheets defining the dune 
% fueld at different times. See example from the North Loup
% River included with this archive. For each dune, we calculate a
% translation (migration) distance in the flow direction, subtract this
% distance and then calculate deformation as the area between the initial
% and final position of the dunes. Other parameters noted in Mason et al. 
% are also measured. 

% Functions finding the distance perpendicular from one bedform to the
% next are adapted from script by Rose Palermo Summer 2015.

% Mackenzie Day
% February 8, 2015
% Updated with additional comments September, 9, 2020

close all; clear all; clc; format long;

% Things to calculate: 1) sinuosity 2) migration speed 3) spur density 4)
% deformation


%% LOAD DATA

% CHANGE THESE PARAMETERS FOR EACH FILE
filename = 'CrestTracingPoints.xls';
num_sheets = 2;%11; %Example data includes 11 sheets
timestamps = [144 149 159 166 173 182 190 200 207 214 219]; % in minutes

% Preallocate variables to hold crest (col 1) and spur (col 2) data
Data = cell(length(num_sheets),2);
Maxes = zeros(num_sheets,3);
Mins = zeros(num_sheets,3);

% Load data from Excel spreadsheet:
for i = 1:num_sheets
    sheet = ['Frame' num2str(i) 'crest'];
    Data(i,1) = {xlsread(filename, sheet, 'a:c')};
    Data(i,2) = {xlsread(filename, sheet, 'e:g')};
    Maxes(i,:) = max(cell2mat(Data(i,1)));
    Mins(i,:) = min(cell2mat(Data(i,1)));
end

% Determine number of crestlines
num_crest = max(Maxes(:,1));
Origx = min(Mins(:,2));
Origy = min(Mins(:,3));

%% ANALYSIS PARAMETERS
search_distance = 50;
search_nodes = 200;
color = 'rbgmck';

%Preallocate space for calculated parameters
DefI = nan(num_crest, num_sheets-1);   %Deformation index
T = nan(num_crest, num_sheets-1);      %Translation distance
PC = nan(num_crest, num_sheets-1);     %Procrustes distance
S = nan(num_crest, num_sheets-1);      %Sinuosity
dS = nan(num_crest, num_sheets-1);     %Change in sinuosity
h1 = figure;
h2 = figure;
h3 = figure;
h4 = figure;
all1 = figure;

%% ITERATE THROUGH FRAMES
for f = 1:num_sheets-1
    F1 = cell2mat(Data(f,1))';
    F2 = cell2mat(Data(f+1,1))';
    dt = timestamps(f+1)-timestamps(f);
    for c = 1:Maxes(f,1)
        I1 = find (F1 == c);
        J1 = I1+1; K1 = I1+2;
        x1 = F1(J1);  y1 = F1(K1);
        sz = length(x1);
        
        I2 = find (F2 == c);
        J2 = I2+1; K2 = I2+2;
        x2 = F2(J2);  y2 = F2(K2);
        
        if isempty(I2) % Error catch for breaks in crestline numbering
            continue
        elseif isempty(I1)
            continue
        end
        
        
        
        %% MANUALLY ENTER UPSTREAM DIRECTION (only once)
        if c == 1 && f ==1
            hc= figure; %make a map that shows the vectors used to calculate the disances
            plot(x1,y1)
            axis equal
            title('Dunes')
            xlabel('Easting')
            ylabel('Northing')
            [xsea,ysea]=ginput(1);
            close(hc);
        end
        
        %Preallocate migration vectors
        mig = zeros(size(x1));
        mig_x = zeros(size(x1));
        mig_y = zeros(size(x1));
        mig_s = zeros(size(x1));
        for i = 1:length(x1) % loop over points in bedform
            
            %x values used to calculate line normal to bedform
            x_i = linspace(x1(i)-search_distance, x1(i)+search_distance,search_nodes);
            
            %find normal line using reciprocal slope (downstream direction)
            y_i = (-1)*(x_i - x1(i)) + y1(i);
            %stream flow defined as NW to SE (slope = -1)
            
            %find interections between normal line and deformed bedform
            [Xi,Yi] = intersections(x_i,y_i,x2,y2);
            
            if isempty(Xi) %if the line does not intersect, set the entry to NaN.
                mig(i) = NaN;
                mig_x(i) = NaN;
                mig_y(i) = NaN;
                mig_s(i) = NaN;
            else % if the line does intersect the shoreline, we will accept the
                %first, closest intersection (there may be more than one)
                srx = (Xi(1)-x1(i));
                sry = (Yi(1)-y1(i));
                mig_x(i) = srx;
                mig_y(i) = sry;
                mig_s(i) = sign(dot([srx sry],[xsea-x1(i),ysea-y1(i)]));%direction
                mig(i) = sqrt((Xi(1)-x1(i)).^2 + (Yi(1)-y1(i)).^2); % Euclidean distance
                
            end
            % Correct for translation
            
        end
        x_bar = mean(mig_x(~isnan(mig_x)));
        y_bar = mean(mig_y(~isnan(mig_y)));
        yt = y1 + mig_y - y_bar;
        xt = x1 + mig_x - x_bar;
        
        % Find the area between the initial and deformed curves
        Px = [x1(1:sz-1,1)';x1(2:sz,1)'];
        Px = [Px;xt(2:sz)';xt(1:sz-1)'];
        Py = [y1(1:sz-1,1)';y1(2:sz,1)'];
        Py = [Py;yt(2:sz)';yt(1:sz-1)'];
        A = polyarea(Px, Py);
        Area = sum(A(~isnan(A)));
        xt_real = xt(~isnan(xt));
        yt_real = yt(~isnan(yt));
        width = sqrt((xt_real(1)-xt_real(end))^2 + (yt_real(1)...
            -yt_real(end))^2);
        
        % Calculate initial, secondary, and change in sinuosity
        curv = [x1 y1];
        [ro co] = size(curv);
        d_calc = pdist2(curv(1:ro-1,:),curv(2:ro,:));
        dis = sum(diag(d_calc));
        w0 = sqrt((x1(1)-x1(end))^2 + (y1(1)-y1(end))^2);
        S1 = dis/w0;
        
        curvt = [xt_real yt_real];
        [rot cot] = size(curvt);
        d_calct = pdist2(curvt(1:rot-1,:),curvt(2:rot,:));
        dis_t = sum(diag(d_calct));
        
        DefI(c,f) = Area/width/dt; %Store deformation index
        T(c,f) = (x_bar^2+y_bar^2)^(1/2); % Store translation distance
        PC(c,f) = (sum(~isnan((xt-x1).^2 +(yt - y1).^2))).^(1/2);
        S(c,f) = dis/width;
        dS(c,f) = S(c,f) - S1;
        
        %% Plot the results
        %make a map that shows the vectors used to calculate the disances
        figure(h1)
        plot(x1-Origx,y1-Origy)
        hold on
        plot(x2-Origx,y2-Origy)
        quiver(x1-Origx,y1-Origy,mig_x,mig_y,0)
        title(['Bedform migration dt = ' ...
            num2str(timestamps(f+1)-timestamps(f)) ' min'])
        xlabel('Easting (m)')
        ylabel('Northing (m)')
        axis equal
        ax = [xlim ylim];
        
        figure(h2)
        hold on
        patch(Px-Origx,Py-Origy,color(1)) %plot polygons
        
        title('Bedform areal changes used to calculate DI')
        xlabel('Easting (m)')
        ylabel('Northing (m)')
        axis equal
        axis(ax)
        
    end
    
    figure(h3)
    scatter(T(:,f), DefI(:,f))
    hold on
    title(['Frame' num2str(f) ' to Frame' num2str(f+1)])
    ylabel('Deformation index (speed?) (m/min)')
    xlabel('Translation distance (m)')
    
    figure(h4)
    scatter(T(:,f), dS(:,f))
    title(['Frame' num2str(f) ' to Frame' num2str(f+1)])
    ylabel('change in sinuosity')
    xlabel('Translation distance (m)')
    
%     saveas(h1, ['Frame' num2str(f) '-' num2str(f+1) '_migration'], 'png')
%     saveas(h2, ['Frame' num2str(f) '-' num2str(f+1) '_polygons'], 'png')
%     saveas(h3, ['T_DefI_Frame' num2str(f) '-' num2str(f+1)], 'png')
%     saveas(h4, ['T_dS_Frame' num2str(f) '-' num2str(f+1)], 'png')
%     close([h1 h2 h3 h4])
    
    figure(all1)
    loglog(T(:,f), DefI(:,f),'.k')
    hold on
    title(['All data Frame' num2str(f) ' to Frame' num2str(f+1)])
    ylabel('Deformation index (speed?) (m/min)')
    xlabel('Translation distance (m)')
end
%saveas(all1, 'AllData_DefI', 'png')
