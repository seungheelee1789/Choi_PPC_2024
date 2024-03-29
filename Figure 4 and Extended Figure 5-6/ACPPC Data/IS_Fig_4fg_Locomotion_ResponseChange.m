% This code calculate Pearson's correlation between locomotion speed and
% individual calcium activities at locomotion onset and offset.
% Run this code in the folder where the exmaple data are.

clear variables; close all; clc;

Mat = FindMatFiles();

MergedOnsetR = []; MergedOnsetP = [];
MergedOffsetR = []; MergedOffsetP = [];

for i = 1:5
    for j = 1:4
       MergedF{i,j} = [];
       MergedAmp{i,j} = [];
       MergedH{i,j} = [];
       
       
       S_MergedF{i,j} = []; S_MergedAmp{i,j} = []; S_MergedH{i,j} = []; S_MergedIndAmp{i,j} = [];
       M_MergedF{i,j} = []; M_MergedAmp{i,j} = []; M_MergedH{i,j} = []; M_MergedIndAmp{i,j} = [];
    end
end

S_MergedInconF{1} = []; S_MergedInconF{2} = [];
M_MergedInconF{1} = []; M_MergedInconF{2} = [];

AmpPost = 1;
AmpPre = 0.2;
PValue = 0.01;

MergedOnsetF = []; MergedOffsetF = [];
MergedOnsetSpeed = []; MergedOffsetSpeed = [];

for MatIdx = 1:numel(Mat)
    load(Mat{MatIdx});
    
    nCell = size(EventF{1,1},2);
    OnsetF = []; OffsetF = [];
    OnsetSpeed = []; OffsetSpeed = [];

    if isempty(LocomotionOffsetF)
    else
        DS_Speed = [];
        if isempty(LocomotionOnsetF) == 0
            DS_Factor = size(LocomotionOnsetSpeed,2)/size(LocomotionOnsetF,1);
            for i = 1:size(LocomotionOnsetF,1)
                DS_Speed(:,i) = nanmean(LocomotionOnsetSpeed(:,(i-1)*DS_Factor+1:i*DS_Factor),2);
            end
            for k = 1:size(LocomotionOnsetF,3)
                OnsetF = cat(1,OnsetF,LocomotionOnsetF(:,:,k));
                OnsetSpeed = cat(1,OnsetSpeed,DS_Speed(k,:)');
            end
            MergedOnsetSpeed = cat(1,MergedOnsetSpeed,nanmean(DS_Speed,1));
            MergedOnsetF = cat(2,MergedOnsetF,nanmean(LocomotionOnsetF,3));
        else
            MergedOnsetSpeed = cat(1,MergedOnsetSpeed,nan(1,ImgHz*4));
            MergedOnsetF = cat(2,MergedOnsetF,nan(ImgHz*4,nCell));
        end

        DS_Speed = [];
        if isempty(LocomotionOffsetF) == 0
            DS_Factor = size(LocomotionOffsetSpeed,2)/size(LocomotionOffsetF,1);
            for i = 1:size(LocomotionOffsetF,1)
                DS_Speed(:,i) = nanmean(LocomotionOffsetSpeed(:,(i-1)*DS_Factor+1:i*DS_Factor),2);
            end
            for k = 1:size(LocomotionOffsetF,3)
                OffsetF = cat(1,OffsetF,LocomotionOffsetF(:,:,k));
                OffsetSpeed = cat(1,OffsetSpeed,DS_Speed(k,:)');
            end
            MergedOffsetSpeed = cat(1,MergedOffsetSpeed,nanmean(DS_Speed,1));
            MergedOffsetF = cat(2,MergedOffsetF,nanmean(LocomotionOffsetF,3));
        else
            MergedOffsetSpeed = cat(1,MergedOffsetSpeed,nan(1,ImgHz*4));
            MergedOffsetF = cat(2,MergedOffsetF,nan(ImgHz*4,nCell));
        end
    end
    for j = 1:nCell
        if isempty(OnsetF)
            MergedOnsetR = cat(1,MergedOnsetR,nan);
            MergedOnsetP = cat(1,MergedOnsetP,nan);
        else
            [R, P] = corrcoef(OnsetF(:,j),OnsetSpeed);
            MergedOnsetR = cat(1,MergedOnsetR,R(1,2));
            MergedOnsetP = cat(1,MergedOnsetP,P(1,2));
        end
        if isempty(OffsetF)
            MergedOffsetR = cat(1,MergedOffsetR,nan);
            MergedOffsetP = cat(1,MergedOffsetP,nan);
        else
            [R, P] = corrcoef(OffsetF(:,j),OffsetSpeed);
            MergedOffsetR = cat(1,MergedOffsetR,R(1,2));
            MergedOffsetP = cat(1,MergedOffsetP,P(1,2));
        end
    end
    

    for i = 1:4
        for j = 1:4
            clearvars H IndAmp
            if isempty(EventF{i,j}) == 0
                AvgF = nanmean(EventF{i,j},3);
                AvgAmp = (nanmean(nanmean(EventF{i,j}(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,1:nCell,:),1)-nanmean(EventF{i,j}((PSTHPre-AmpPre)*ImgHz+1:PSTHPre*ImgHz,1:nCell,:),1),3));
                for N = 1:nCell
                    clearvars p
                    IndAmp{i,j}{N}(:,1) = nanmean(EventF{i,j}(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,N,:),1)-nanmean(EventF{i,j}((PSTHPre-AmpPre)*ImgHz+1:PSTHPre*ImgHz,N,:),1);
                    p = signrank(IndAmp{i,j}{N}(:,1));
                    if p > PValue
                        H(1,N) = 0;
                    elseif p < PValue & mean(IndAmp{i,j}{N}(:,1)) > 0
                        H(1,N) = 1;
                    elseif p < PValue & mean(IndAmp{i,j}{N}(:,1)) < 0
                        H(1,N) = -1;
                    end
                end
            else

                clearvars AvgF AvgAmp
                AvgF(1:(PSTHPre+PSTHPost)*ImgHz,1:nCell) = nan;
                AvgAmp(1,1:nCell) = nan;    
                H(1,1:nCell) = nan;
                for N = 1:nCell
                    IndAmp{i,j}{N} = [];
                end
            end
            MergedF{i,j} = cat(2,MergedF{i,j},AvgF);
            MergedAmp{i,j} = cat(2,MergedAmp{i,j},AvgAmp);
            MergedH{i,j} = cat(2,MergedH{i,j},H);
            if strcmp(Mat{MatIdx}(1:6),'Moving') == 1
                M_MergedF{i,j} = cat(2,M_MergedF{i,j},AvgF);
                M_MergedAmp{i,j} = cat(2,M_MergedAmp{i,j},AvgAmp);
                M_MergedH{i,j} = cat(2,M_MergedH{i,j},H);
                M_MergedIndAmp{i,j} = cat(2,M_MergedIndAmp{i,j},IndAmp{i,j});
            elseif strcmp(Mat{MatIdx}(1:10),'Stationary') == 1
                S_MergedF{i,j} = cat(2,S_MergedF{i,j},AvgF);
                S_MergedAmp{i,j} = cat(2,S_MergedAmp{i,j},AvgAmp);
                S_MergedH{i,j} = cat(2,S_MergedH{i,j},H);
                S_MergedIndAmp{i,j} = cat(2,S_MergedIndAmp{i,j},IndAmp{i,j});
            else
                error('There is something wrong');
            end
        end
        if i == 4
            CatData = cat(3,EventF{i,1},EventF{i,2});
            AvgAmp = (nanmean(nanmean(CatData(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,1:nCell,:),1)-nanmean(CatData((PSTHPre-AmpPre)*ImgHz+1:PSTHPre*ImgHz,1:nCell,:),1),3));
            InconF = nanmean(CatData,3);
            if strcmp(Mat{MatIdx}(1:10),'Stationary') == 1
                S_MergedF{5,1} = cat(2,S_MergedF{5,1},InconF);
                S_MergedAmp{5,1} = cat(2,S_MergedAmp{5,1},AvgAmp);
            elseif strcmp(Mat{MatIdx}(1:6),'Moving') == 1
                M_MergedF{5,1} = cat(2,M_MergedF{5,1},InconF);
                M_MergedAmp{5,1} = cat(2,M_MergedAmp{5,1},AvgAmp);
            end
            CatData = cat(3,EventF{i,3},EventF{i,4});
            AvgAmp = (nanmean(nanmean(CatData(PSTHPre*ImgHz+1:(PSTHPre+AmpPost)*ImgHz,1:nCell,:),1)-nanmean(CatData((PSTHPre-AmpPre)*ImgHz+1:PSTHPre*ImgHz,1:nCell,:),1),3));
            InconF = nanmean(CatData,3);
            if strcmp(Mat{MatIdx}(1:10),'Stationary') == 1
                S_MergedF{5,3} = cat(2,S_MergedF{5,3},InconF);
                S_MergedAmp{5,3} = cat(2,S_MergedAmp{5,1},AvgAmp);
            elseif strcmp(Mat{MatIdx}(1:6),'Moving') == 1
                M_MergedF{5,3} = cat(2,M_MergedF{5,3},InconF);
                M_MergedAmp{5,3} = cat(2,M_MergedAmp{5,3},AvgAmp);
            end
        end
    end
end

%
PValue2 = 0.01;

Target = [1 1; 1 3; 2 1; 2 3];

for i = 1:size(Target,1)
    for ii = 1:5
        Change{ii,i} = []; % 1 Same; 2 Stronger; 3 Weaker; 4 Reversed; 5 No response
    end

    x = Target(i,1); y= Target(i,2);

    Change{5,i} = find(S_MergedH{x,y} == 0 & M_MergedH{x,y} == 0);
    Idx{1,i} = find(S_MergedH{x,y} ~= 0 | M_MergedH{x,y} ~= 0);
    S_Amp = S_MergedIndAmp{x,y}(Idx{i});
    M_Amp = M_MergedIndAmp{x,y}(Idx{i});

    for j = 1:numel(Idx{i});
        if S_MergedH{x,y}(Idx{i}(j)) == 0 & M_MergedH{x,y}(Idx{i}(j)) ~= 0
            Change{2,i} = cat(2,Change{2,i},Idx{i}(j));
        elseif S_MergedH{x,y}(Idx{i}(j)) ~= 0 & M_MergedH{x,y}(Idx{i}(j)) == 0
            Change{3,i} = cat(2,Change{3,i},Idx{i}(j));
        elseif S_MergedH{x,y}(Idx{i}(j)).*M_MergedH{x,y}(Idx{i}(j)) == -1
            Change{4,i} = cat(2,Change{4,i},Idx{i}(j));
        elseif S_MergedH{x,y}(Idx{i}(j)) ~= 0 & M_MergedH{x,y}(Idx{i}(j)) ~= 0
            p = ranksum(S_Amp{j},M_Amp{j});
            %[~,p] = ttest2(S_Amp{j},M_Amp{j});
            if p > PValue2
                Change{1,i} = cat(2,Change{1,i},Idx{i}(j));
            elseif p < PValue2
                if abs(nanmean(S_Amp{j})) > abs(nanmean(M_Amp{j}))
                    Change{3,i} = cat(2,Change{3,i},Idx{i}(j));
                else
                    Change{2,i} = cat(2,Change{2,i},Idx{i}(j));
                end
            end
        end
    end
end

%

mkdir('Figure');
cd('Figure');


xAxis = -2+1/ImgHz:1/ImgHz:2;

FAxis = [-2 2];
AllColor = [0 180 0]./255;

close all;

fig = figure('Position', [50 50 90 110]); hold on;

AllF = MergedOnsetF;
AllF_Avg = nanmean(AllF,2); AllF_Sem = nanstd(AllF,0,2)./sqrt(size(AllF,2)-numnan(AllF(1,:)));
Speed_Avg = nanmean(MergedOnsetSpeed,1); Speed_Sem = nanstd(MergedOnsetSpeed,1)./sqrt(size(MergedOnsetSpeed,1)-numnan(MergedOnsetSpeed(:,1)));

xlim([-1 2]);

yyaxis right
errorshade(xAxis,Speed_Avg,Speed_Avg+Speed_Sem,Speed_Avg-Speed_Sem,[0 0 0]);
alpha(0.2);
plot(xAxis,Speed_Avg,'lineWidth',1,'color','k');
ylabel('Speed (cm/s)','FontName','Arial','FontSize',6);
xlim([-1 2]);

set(gca,'YColor',[0 0 0]);

yyaxis left

errorshade(xAxis,AllF_Avg',AllF_Avg'+AllF_Sem',AllF_Avg'-AllF_Sem',AllColor);
alpha(0.2);

plot(xAxis,AllF_Avg,'lineWidth',1,'color',AllColor,'lineStyle','-');

set(gca,'YColor',[0 180 0]./255);

%xlim([-1 2]); %ylim([-1 1]); %yticks([-1 0 1]);
ylim([-1 1]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Time (s)','FontName','Arial','FontSize',6);
ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

saveas(fig,['Fig 4f, Locomotion Onset PSTH.svg']);


fig = figure('Position', [50 50 90 110]); hold on;

AllF = MergedOffsetF;
AllF_Avg = nanmean(AllF,2); AllF_Sem = nanstd(AllF,0,2)./sqrt(size(AllF,2)-numnan(AllF(1,:)));
Speed_Avg = nanmean(MergedOffsetSpeed,1); Speed_Sem = nanstd(MergedOffsetSpeed,1)./sqrt(size(MergedOffsetSpeed,1)-numnan(MergedOffsetSpeed(:,1)));

xlim([-1 2]);

yyaxis right
errorshade(xAxis,Speed_Avg,Speed_Avg+Speed_Sem,Speed_Avg-Speed_Sem,[0 0 0]);
alpha(0.2);
plot(xAxis,Speed_Avg,'lineWidth',1,'color','k');
ylabel('Speed (cm/s)','FontName','Arial','FontSize',6);
xlim([-1 2]);

set(gca,'YColor',[0 0 0]);

yyaxis left

errorshade(xAxis,AllF_Avg',AllF_Avg'+AllF_Sem',AllF_Avg'-AllF_Sem',AllColor);
alpha(0.2);

plot(xAxis,AllF_Avg,'lineWidth',1,'color',AllColor,'lineStyle','-');

set(gca,'YColor',[0 180 0]./255);

%xlim([-1 2]); %ylim([-1 1]); %yticks([-1 0 1]);
ylim([-1 1]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Time (s)','FontName','Arial','FontSize',6);
ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

saveas(fig,['Fig 4g, Locomotion Offset PSTH.svg']);

Bin = -1:0.1:1;

fig = figure('Position', [50 50 90 110]); hold on;

%IncR = MergedR(Change{2,i});
AllR = MergedOnsetR;

histogram(AllR,Bin,'FaceColor','none','EdgeColor','k');
scatter(nanmean(AllR),30,'v','MarkerFaceColor','r','MarkerEdgeColor','none');
line([0 0],[0 30],'color','k','lineWidth',1,'lineStyle',':');

xlim([-1 1]);
%ylim(YAxis);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Corr. coef, r','FontName','Arial','FontSize',6);
ylabel('Cell number','FontName','Arial','FontSize',6);

saveas(fig,'Fig 4f, Locomotion Onset Correlation Histogram.svg');

fig = figure('Position', [50 50 90 110]); hold on;

AllR = MergedOffsetR;

histogram(AllR,Bin,'FaceColor','none','EdgeColor','k');
scatter(nanmean(AllR),30,'v','MarkerFaceColor','r','MarkerEdgeColor','none');
line([0 0],[0 30],'color','k','lineWidth',1,'lineStyle',':');

xlim([-1 1]);
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Corr. coef, r','FontName','Arial','FontSize',6);
ylabel('Cell number','FontName','Arial','FontSize',6);

saveas(fig,'Fig 4g, Locomotion Offset Correlation Histogram.svg');

cd ../
