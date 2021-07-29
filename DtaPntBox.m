function DtaPntBox(DataBoxes, Labels, LabelsScatter, ClrSheme, ClrShemeSct, XTickName, yAxisName, Wdth, MkrSiz, RepMean, Lgd, Textex, Save)
%
% Boxplot with possibility of displaying color gradient of
% numerical values of data points together with subgroups coloring.
% 
% Data from example file:
% Numerical value of object response;
% Seven different treatment (A, B, C, D, E, F, G);
% Three different day (1,2,3);
% Twenty one observations overall, seven for subgroup with seven different 
% dosages (1,2,3,4,5,6,7);
%  
% Example boxplot can represent:
% If dosage of compound A in subgroup 2 influences value of response.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT E.G:
%
% DtaPntBox(X, Y, Z, [], Colors, XTickName, 'Rel. int', 1, 32,'Yes','Data points color gradient','yes');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%
% 1) DataBoxes      - Data from e.g. compound concentration with different methods in columns e.g.: Treatment A,
%                     treatment B, (...) treatment F
%
% 2) Lables         - Column vector with Labels ex: Day1, Day2, Day3 - 
%                     This vector indicates how data will be split. 
%                     If you have 3 unique labels it will be split into 3 boxes             
%                     in single Xlabel
%
% 3) LabelsScatter  - Column vector with group annotation for gradient in data
%                     points
%
% 4) ClrSheme       - Color matrix for boxplot filling - it should have the same 
%                     number as unique labels in Labels input.
%               
% 5) ClrShemeSct    - Struct contains cells with matrices with desired colors for data points
%                     representation. This struct should have same number of
%                     cell as unique labels in Labels and each matrix should contain as many
%                     vectors with colors as many unique labels in LabelsScatter
%
% 6) XTickName      - Names for different columns from DataBoxes, e.g.: Treatment A,
%                     treatment B, (...) treatment F
%
% 7) yAxisName      - Label for Y axis
%
% 8) Wdth           - Line width for line plot
%
% 9) MkrSiz         - Marker size for scatter plot. Empty if without markers
%
% 10) RepMean       - Input 'Yes' or 1 for mean group representation
%
% 11) Lgd           - Input 'Yes' or 1 for displaying legend
%
% 12) Textex        - Title of boxplot and name of file if saved
%
% 13) Save          - Input 'Yes' or 1, if you want save graphical output into *.TIF file
%
% Author: Wojciech Wojtowicz 
% Contact: wojciech.wojtowicz@pwr.edu.pl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checking input values sizes

ChckX = size(DataBoxes);
ChckZ = size(Labels);
ChckY = size(LabelsScatter);

if ChckX(1,1) == ChckZ(1,2)
else
    warning(['Check dimension of input data. DataBoxes dimension ' num2str(ChckX(1,1)) 'x' num2str(ChckX(1,2)) ' when Lables '  num2str(ChckZ(1,1)) 'x' num2str(ChckZ(1,2))]) 
return;
end

if ChckX(1,1) == ChckY(1,1) && ChckX(1,2) == ChckY(1,2)
else
    warning(['Check dimension of input data. DataBoxes dimension ' num2str(ChckX(1,1)) 'x' num2str(ChckX(1,2)) ' when LabelsScatter '  num2str(ChckY(1,1)) 'x' num2str(ChckY(1,2))]) 
return;
end
%% Checking RepMean type
if isempty(RepMean)
        RepMean = 0;
elseif isstring(RepMean) || iscell(RepMean) || ischar(RepMean)
    if strcmp(RepMean, 'yes') || strcmp(RepMean, 'Yes')
        RepMean = 1;
    else
        RepMean = 0;
    end
elseif isnumeric(RepMean)
    if RepMean ==1
        RepMean = 1;
    else
        RepMean = 0;
    end
end
 
%% Checking Lgd type
if isempty(Lgd)
        Lgd = 0;
elseif isstring(Lgd) || iscell(Lgd) || ischar(Lgd)
    if strcmp(Lgd, 'yes') || strcmp(Lgd, 'Yes')
        Lgd = 1;
    else
        Lgd = 0;
    end
elseif isnumeric(Lgd)
    if Lgd ==1
        Lgd = 1;
    else
        Lgd = 0;
    end
end
 
%%
[x,~,z] = unique(Labels);
iY =1;

for i=1:length(x)    
    Nums(i) = length(x)+1-i;    
end

figure('units','normalized','outerposition',[0 0 1 0.7]); set(gcf,'color','w');
box on; grid on;

for k=1:size(DataBoxes,2)
    Spt = 0;
    j= 1;
for i=1:length(Nums)
 
Indxs = find(z==x(i)); 
hold on
%% Data points ploting
[x1,~,~] = unique(LabelsScatter(:,i));
if isempty(ClrShemeSct)
    ClrSct = repmat([0], length(x),length(z)/length(x))';
else
    ClrSct = cell2mat(ClrShemeSct(i));
end


for h=1:length(x1)
SctDta = DataBoxes(Indxs,k);
IndxDit = find(LabelsScatter(Indxs)==x1(h));
if isempty(MkrSiz)
else
scatter(repmat(mean([k+Spt, k+Spt+0.25]),1, length(SctDta(IndxDit,1))), SctDta(IndxDit,1), MkrSiz, ClrSct(h,:), 'filled','MarkerFaceAlpha',0.85','jitter','on','jitterAmount',0.05); %'o','MarkerEdgeColor',ClrSheme(i,:), 'MarkerFaceColor', ClrSheme(i,:))
end
if h==1
if RepMean == 1  
   MeanForPlt(iY) = mean(SctDta);
   MeanForPltY(iY) = mean([k+Spt, k+Spt+0.25]);
   iY = iY+1;
end
end  

end

%% Box
plot([k+Spt, k+Spt+0.25], [quantile(DataBoxes(Indxs,k),0.25), quantile(DataBoxes(Indxs,k),0.25)],'k','LineWidth', Wdth)
plot([k+Spt, k+Spt+0.25], [quantile(DataBoxes(Indxs,k),0.75), quantile(DataBoxes(Indxs,k),0.75)],'k','LineWidth', Wdth)
plot([k+Spt, k+Spt],[quantile(DataBoxes(Indxs,k),0.25), quantile(DataBoxes(Indxs,k),0.75)],'k','LineWidth', Wdth)
plot([k+Spt+0.25, k+Spt+0.25],  [quantile(DataBoxes(Indxs,k),0.25), quantile(DataBoxes(Indxs,k),0.75)],'k','LineWidth', Wdth)
 
%% Whisker
plot([mean([k+Spt, k+Spt+0.25]), mean([k+Spt, k+Spt+0.25])], [quantile(DataBoxes(Indxs,k),0.25), min(DataBoxes(Indxs,k))],'--k','LineWidth', Wdth)
plot([mean([k+Spt, k+Spt+0.25]), mean([k+Spt, k+Spt+0.25])], [quantile(DataBoxes(Indxs,k),0.75), max(DataBoxes(Indxs,k))],'--k','LineWidth', Wdth)
plot([mean([k+Spt, k+Spt+0.25])-0.05, mean([k+Spt, k+Spt+0.05])+0.15], [min(DataBoxes(Indxs,k)),min(DataBoxes(Indxs,k))],'k','LineWidth', Wdth)
plot([mean([k+Spt, k+Spt+0.25])-0.05, mean([k+Spt, k+Spt+0.05])+0.15], [max(DataBoxes(Indxs,k)),max(DataBoxes(Indxs,k))],'k','LineWidth', Wdth)

%% Mean bar
plot([k+Spt, k+Spt+0.25], [mean(DataBoxes(Indxs,k)),mean(DataBoxes(Indxs,k))],'k','LineWidth', Wdth)
 
%% Box coloring
if isempty(ClrSheme)
else
a1 = [k+Spt,quantile(DataBoxes(Indxs,k),0.25); k+Spt+0.25, quantile(DataBoxes(Indxs,k),0.25); k+Spt+0.25, quantile(DataBoxes(Indxs,k),0.75); k+Spt, quantile(DataBoxes(Indxs,k),0.75)];
f = [1 2 3 4];
patch('Faces',f,'Vertices',a1,'FaceColor',ClrSheme(i,:),'FaceAlpha',.25);
end
Spt = Spt+0.25;
clear Indxs

end
end

%% Setting labels for XTick
set(gca, 'XTick', [1+length(x)*0.25/2:1:size(DataBoxes,2)+0.25*length(Nums)]);
xlim([min(x)-0.125, size(DataBoxes,2)+0.25*max(x)+0.125]);
xticklabels(XTickName);


%% Ploting vertical lines between boxes
NoVlines = size(XTickName,2);
Num = NoVlines;
IdYim = ylim;

for i=1:NoVlines
if Num > 1
NumForVline =  Num + ((size(DataBoxes,2)+0.25*length(x)) - (size(DataBoxes,2)+0.25*length(x)+0.25))/2;
plot([NumForVline, NumForVline], IdYim, 'k', 'LineWidth',0.25)
Num =  Num - 1;
end
end
SpltY=1;
if RepMean == 1 
for i = 1:size(DataBoxes,2)
    plot(MeanForPltY(SpltY:SpltY+(length(x)-1)), MeanForPlt(SpltY:SpltY+(length(x)-1)),'--ok', 'LineWidth', Wdth,'MarkerFaceColor','k');
    SpltY = SpltY+length(x);
end
end

%% Visiual aspect
ylim([IdYim(1), IdYim(2)]);
ylabel(yAxisName);
set(gca, 'FontSize',20, 'FontWeight','Bold');
%% Set title
if isempty(Textex)
else
title(Textex)
end

%% Set legend
if Lgd == 1
    if isempty(ClrSheme)
    else
for i=1:length(x)
    if isnumeric(x(i))
        xIn = num2str(x(i));
    else
        xIn = x(i);
    end
%          CmbFig(i) = plot(NaN,NaN,'s','Color',ClrSheme(i,:),'MarkerFaceColor',ClrSheme(i,:),'MarkerSize',16,'DisplayName', xIn);
%          CmbFig(i) = plot(NaN,NaN,'s','Color',ClrSheme(i,:),'LineWidth',9,'DisplayName', xIn);
%          CmbFig(i).Color(4)=0.25;
   CmbFig(i) = scatter(NaN,NaN,10,'s','LineWidth',10,'MarkerFaceColor',ClrSheme(i,:),'MarkerEdgeColor',ClrSheme(i,:),'MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25,'DisplayName', xIn);
%          CmbFig(i).SizeData = 60;
end
        legend(CmbFig,'Location','eastoutside');
        
%  [~, MkrData] = legend(CmbFig,'Location','eastoutside');
%  for i=1:length(x)
%  ObjMkrData = findobj(MkrData(i+length(x)), 'type', 'patch'); % objects of legend of type patch
%  ObjMkrData = findobj(MkrData(i+length(x)), 'type','line')
%  set(ObjMkrData, 'Markersize', 25); % set marker size as desired
%  set(ObjMkrData, 'FaceAlpha', 0.25);
    end
end
%% Check Save input style
if isempty(Save)
        Save = 0;
elseif isstring(Save) || iscell(Save) || ischar(Save)
    if strcmp(Save, 'yes') || strcmp(Save, 'Yes')
        Save = 1;
    else
        Save = 0;
    end
elseif isnumeric(Save)
    if Save ==1
        Save = 1;
    else
        Save = 0;
    end
end
if Save == 1
if isempty(Textex)
    Textex='';
end
   print([ 'Boxplot-' Textex],'-dtiff','-r0');
end
end