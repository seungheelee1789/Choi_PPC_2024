% You should have run "IS_Run_dPCA.m" before this
% This code makes figures about projection distance and Go/Nogo distance
% Run this code in the folder that has "dPCA Result.mat" file

clear all; close all; clc;

SubSpace = 1; % Subspace type, 1 = 'Stimulus'; 2 = 'Decision'; 3 = 'Condition-independent'; 4 = 'Stim-Deci Interaction';

UniTarget = [1 4 5 8];
InconTarget = [2 3; 2 3; 1 4; 1 4];
ConTarget = [1 3; 2 4];
V_Target{3,1} = [1 1]; V_Target{3,2} = [1 2]; V_Target{4,1} = [1 2]; V_Target{4,2} = [1 1];
A_Target{3,1} = [2 1]; A_Target{3,2} = [2 2]; A_Target{4,1} = [2 1]; A_Target{4,2} = [2 2];

DistTime1 = 0; 
DistTime2 = 1; % Time point from the stimulus onset  to calculate Euclidean distance (in second)
%DistBin = (Pre+DistTime1)*ImgHz:(Pre+DistTime2)*ImgHz;

PC = 1:3 %PCs to use for distance calculation, remove this line to use all PCs (for Fig S3g)

load('dPCA Result.mat');
N = 0; nS = 0; nM = 0; SIdx = []; MIdx = [];
for Idx = 1:size(AvgScore,1)
    if isempty(AvgScore{Idx,SubSpace}) == 0
        N = N+1;
        State{N} = Name{Idx,2};
        for j = 1:numel(AvgScore{Idx,SubSpace})
            ScoreData{N}{j} = [];
            for jj = 1:size(AvgScore{Idx,SubSpace}{j},2)
                ScoreData{N}{j} = cat(2,ScoreData{N}{j},AvgScore{Idx,SubSpace}{j}(:,jj));
            end
        end
        if strcmp(Name{Idx,2},'m')
            for i = 1:4
                for j = 1:2
                    for k = 1:2
                        MvmtData{N,1}{i,j}{k,1} = MvmtSortedScore{Idx,SubSpace}{i,j}{k};
                    end
                    MvmtData{N,1}{i,j}{3,1} = [];
                end
                for k = 1:2
                    InconMvmtData{N,1}{1,i}{k,1} = InconMvmtSortedScore{Idx,SubSpace}{i}{k};
                end
                InconMvmtData{N,1}{1,i}{3,1} = [];
            end
        else
            for i = 1:4
                for j = 1:2
                    for k = 1:2
                        MvmtData{N,1}{i,j}{k,1} = [];
                    end
                    k = 3;
                    MvmtData{N,1}{i,j}{3,1} = MvmtSortedScore{Idx,SubSpace}{i,j}{k};
                end
                k = 3;
                InconMvmtData{N,1}{1,i}{3,1} = InconMvmtSortedScore{Idx,SubSpace}{i}{k};
            end
        end
    end
    if strcmp(Name{Idx,2},'m'); nM = nM + 1; MIdx = cat(1,MIdx,Idx);
    elseif strcmp(Name{Idx,2},'s'); nS = nS + 1; SIdx = cat(1,SIdx,Idx); end;
end

M_Session = [];
S_Session = [];
for Idx = 1:numel(State)
    if strcmp(State{Idx},'m') == 1
        M_Session = cat(1,M_Session,Idx);
        for i = 1:2
            ADomRate(1,i,Idx) = InconTrialNumber{Idx}(i,1)/sum(InconTrialNumber{Idx}(i,1:2));
            ADomRate(2,i,Idx) = InconTrialNumber{Idx}(i,3)/sum(InconTrialNumber{Idx}(i,3:4));
        end
    elseif strcmp(State{Idx},'s') == 1
        S_Session = cat(1,S_Session,Idx);
        ADomRate(1,3,Idx-nM) = InconTrialNumber{Idx}(3,1)/sum(InconTrialNumber{Idx}(3,1:2));
        ADomRate(2,3,Idx-nM) = InconTrialNumber{Idx}(3,3)/sum(InconTrialNumber{Idx}(3,3:4));
    end
end

Incon_Dv = []; Incon_Da = [];
Con_Dv = []; Con_Da = [];
Uni_D = [];
for Idx = 1:numel(ScoreData)
    CatUni = []; RMS = 0;
    for i = 1:numel(UniTarget)
        U_i = UniTarget(i);
        UniData{i} = ScoreData{Idx}{U_i};
        CatUni = cat(1,CatUni,UniData{i});
    end

    if exist('PC')
        if size(UniData{1},2) < max(PC)
            PC = 1:size(UniData{1},2);
        end
    else
        disp('There is no PC input!!!!!!!!!!');
        PC = 1:size(UniData{1},2);
    end

    %RMS = rms(CatUni(:,PC),'all');
    RMS = 1;

    for i = 13:16
        InconData{i-12} = ScoreData{Idx}{i};
    end

    for i = 1:numel(InconData)
        V_i = InconTarget(i,1); A_i = InconTarget(i,2);
        for ii = 1:(Pre+Post)*ImgHz+1
            Incon_Dv(ii,i,Idx) = (norm(InconData{i}(ii,PC) - UniData{V_i}(ii,PC)))./RMS;
            Incon_Da(ii,i,Idx) = (norm(InconData{i}(ii,PC) - UniData{A_i}(ii,PC)))./RMS;
        end
        InconRatio(:,i,Idx) = (Incon_Da(:,i,Idx) - Incon_Dv(:,i,Idx))./(Incon_Da(:,i,Idx) + Incon_Dv(:,i,Idx));
    end

    ConData{1} = ScoreData{Idx}{9};
    ConData{2} = ScoreData{Idx}{12};
    for i = 1:numel(ConData)
        V_i = ConTarget(i,1); A_i = ConTarget(i,2);
        for ii = 1:(Pre+Post)*ImgHz+1
            Con_Dv(ii,i,Idx) = norm(ConData{i}(ii,PC) - UniData{V_i}(ii,PC))./RMS;
            Con_Da(ii,i,Idx) = norm(ConData{i}(ii,PC) - UniData{A_i}(ii,PC))./RMS;
        end
        ConRatio(:,i,Idx) = (Con_Da(:,i,Idx) - Con_Dv(:,i,Idx))./(Con_Da(:,i,Idx) + Con_Dv(:,i,Idx));
    end

    for x = 3:4
        for y = 1:2
            V_x = V_Target{x,y}(1); V_y = V_Target{x,y}(2);
            A_x = A_Target{x,y}(1); A_y = A_Target{x,y}(2);
            for z = 1:3
                if isempty(MvmtData{Idx}{x,y}{z})
                else
                    for ii = 1:(Pre+Post)*ImgHz+1
                        if z < 3
                            MvmtDv{x,y}{z}(ii,Idx) = norm(MvmtData{Idx}{x,y}{z}(ii,PC) - MvmtData{Idx}{V_x,V_y}{z}(ii,PC))./RMS;
                            MvmtDa{x,y}{z}(ii,Idx) = norm(MvmtData{Idx}{x,y}{z}(ii,PC) - MvmtData{Idx}{A_x,A_y}{z}(ii,PC))./RMS;
                        else
                            MvmtDv{x,y}{z}(ii,Idx-nM) = norm(MvmtData{Idx}{x,y}{z}(ii,PC) - MvmtData{Idx}{V_x,V_y}{z}(ii,PC))./RMS;
                            MvmtDa{x,y}{z}(ii,Idx-nM) = norm(MvmtData{Idx}{x,y}{z}(ii,PC) - MvmtData{Idx}{A_x,A_y}{z}(ii,PC))./RMS;
                        end
                    end
                end
            end
        end
    end

    for x = 1:4
        if x < 3
            V_i = 2; A_i = 1;
        else
            V_i = 1;A_i = 2;
        end

        for z = 1:3
            if isempty(InconMvmtData{Idx}{x}{z}) == 0
                if z < 3
                    for ii = 1:(Pre+Post)*ImgHz+1
                        MvmtInconDv{x,z}(ii,Idx) = norm(InconMvmtData{Idx}{x}{z}(ii,PC) - MvmtData{Idx}{1,V_i}{z}(ii,PC))./RMS;
                        MvmtInconDa{x,z}(ii,Idx) = norm(InconMvmtData{Idx}{x}{z}(ii,PC) - MvmtData{Idx}{2,A_i}{z}(ii,PC))./RMS;
                    end
                else
                    for ii = 1:(Pre+Post)*ImgHz+1
                        MvmtInconDv{x,z}(ii,Idx-nM) = norm(InconMvmtData{Idx}{x}{z}(ii,PC) - MvmtData{Idx}{1,V_i}{z}(ii,PC))./RMS;
                        MvmtInconDa{x,z}(ii,Idx-nM) = norm(InconMvmtData{Idx}{x}{z}(ii,PC) - MvmtData{Idx}{2,A_i}{z}(ii,PC))./RMS;
                    end
                end
            end
        end
    end

    for ii = 1:(Pre+Post)*ImgHz+1
        Uni_D(ii,1,Idx) = (norm(UniData{1}(ii,PC) - UniData{2}(ii,PC)))/RMS;
        Uni_D(ii,2,Idx) = (norm(UniData{1}(ii,PC) - UniData{3}(ii,PC)))/RMS;
        Uni_D(ii,3,Idx) = (norm(UniData{1}(ii,PC) - UniData{4}(ii,PC)))/RMS;
        Uni_D(ii,4,Idx) = (norm(UniData{2}(ii,PC) - UniData{3}(ii,PC)))/RMS;
        Uni_D(ii,5,Idx) = (norm(UniData{2}(ii,PC) - UniData{4}(ii,PC)))/RMS;
        Uni_D(ii,6,Idx) = (norm(UniData{3}(ii,PC) - UniData{4}(ii,PC)))/RMS;
        Uni_D(ii,7,Idx) = (norm(ConData{1}(ii,PC) - ConData{2}(ii,PC)))/RMS;
        Uni_D(ii,8,Idx) = (norm(InconData{3}(ii,PC) - InconData{2}(ii,PC)))/RMS;
        Uni_D(ii,9,Idx) = (norm(InconData{4}(ii,PC) - InconData{1}(ii,PC)))/RMS;
    end

    for x = 1:3
        for y = 1:2
            for z = 1:3
                if isempty(MvmtData{Idx}{x,y}{z}) == 0
                    if z < 3;
                        PD{x,y}(z,Idx) = norm(MvmtData{Idx}{x,y}{z}((Pre+DistTime2)*ImgHz,PC) - MvmtData{Idx}{x,y}{z}(Pre*ImgHz,PC))./RMS;
                    else
                        PD{x,y}(z,Idx-nM) = norm(MvmtData{Idx}{x,y}{z}((Pre+DistTime2)*ImgHz,PC) - MvmtData{Idx}{x,y}{z}(Pre*ImgHz,PC))./RMS;
                    end
                end
            end
        end
        for z = 1:3
            if isempty(MvmtData{Idx}{x,y}{z}) == 0
                if z < 3;
                    for ii = 1:(Pre+Post)*ImgHz+1
                        GNG_D{x,1}{z}(ii,Idx) = norm(MvmtData{Idx}{x,1}{z}(ii,PC) - MvmtData{Idx}{x,2}{z}(ii,PC))./RMS;
                    end
                else
                    for ii = 1:(Pre+Post)*ImgHz+1
                        GNG_D{x,1}{z}(ii,Idx-nM) = norm(MvmtData{Idx}{x,1}{z}(ii,PC) - MvmtData{Idx}{x,2}{z}(ii,PC))./RMS;
                    end
                end
            end
        end
    end
    for z = 1:3
        if isempty(MvmtData{Idx}{4,1}{z}) == 0
            if z < 3
                for ii = 1:(Pre+Post)*ImgHz+1
                    GNG_D{4,1}{z}(ii,Idx) = norm(InconMvmtData{Idx}{1,1}{z}(ii,PC) - InconMvmtData{Idx}{1,3}{z}(ii,PC))./RMS;
                    GNG_D{4,2}{z}(ii,Idx) = norm(InconMvmtData{Idx}{1,2}{z}(ii,PC) - InconMvmtData{Idx}{1,4}{z}(ii,PC))./RMS;
                end
            else
                for ii = 1:(Pre+Post)*ImgHz+1
                    GNG_D{4,1}{z}(ii,Idx-nM) = norm(InconMvmtData{Idx}{1,1}{z}(ii,PC) - InconMvmtData{Idx}{1,3}{z}(ii,PC))./RMS;
                    GNG_D{4,2}{z}(ii,Idx-nM) = norm(InconMvmtData{Idx}{1,2}{z}(ii,PC) - InconMvmtData{Idx}{1,4}{z}(ii,PC))./RMS;
                end
            end
        end
    end

    for x = 1:4
        for z = 1:3
            if isempty(InconMvmtData{Idx}{x}{z}) == 0
                if z < 3
                    PD{4,x}(z,Idx) = norm(InconMvmtData{Idx}{x}{z}((Pre+DistTime2)*ImgHz,PC) - InconMvmtData{Idx}{x}{z}(Pre*ImgHz,PC))./RMS;
                else
                    PD{4,x}(z,Idx-nM) = norm(InconMvmtData{Idx}{x}{z}((Pre+DistTime2)*ImgHz,PC) - InconMvmtData{Idx}{x}{z}(Pre*ImgHz,PC))./RMS;
                end
            end
        end
    end
end

%%
PlotData = [];
for Idx = 1:numel(ComponentInfo)
    Data = ComponentInfo{Idx};
    for i = 1:4
        PlotData{i}(1:20,Idx) = nan;
        PlotData{i}(1:numel(find(Data.whichMarg == i)),Idx) = (Data.componentVar(find(Data.whichMarg == i)))';
    end
end

Color = [207 113 174; 193 231 227; 180 180 180; 174 207 113]./255;

for State = 1:3
    if State == 1; Session = SIdx; elseif State == 2; Session = MIdx; elseif State == 3; Session = [MIdx; SIdx;]; end;
    
    Fig = figure('Position',[50 50 82 110]);
    hold on
    for i = 1:4
        Avg = []; Sem = [];
        Avg = nanmean(PlotData{i}(:,Session),2);
        for ii = 1:size(PlotData{1},1)
            Sem(ii,1) = nanstd(PlotData{i}(ii,Session),0,2)./sqrt(numel(Session)-numnan(PlotData{i}(ii,Session)));
        end
        xTicks = 1:(20-numnan(Avg));
        End = 20-numnan(Avg);
        errorshade(xTicks,Avg(1:End)',Avg(1:End)'+Sem(1:End)',Avg(1:End)'-Sem(1:End)',Color(i,:));
        alpha(0.5);
    end
    for i = 1:4
        Avg = [];
        Avg = nanmean(PlotData{i}(:,Session),2);
        plot(Avg,'lineWidth',1,'color',Color(i,:));
    end
    xlim([0.5 8.5]); xticks([1:9]);
    ylim([0 40]);
    %ylim([0 5]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('PC#','FontName','Arial','FontSize',6);
    ylabel('Explained variance (%)','FontName','Arial','FontSize',6);
    
    mkdir('Figure'); cd('Figure');
    if State == 1;
        saveas(Fig,'Fig S2e, Stationary Explained Variance.svg');
    elseif State == 2
        saveas(Fig,'Fig S2e, Moving Explained Variance.svg');
    elseif State == 3
        saveas(Fig,'Fig 2f, Total Explained Variance.svg');
    end
    cd ../
end
%% ANOVA test
% AnovaIdx = 1:10; AnovaPCIdx = 4;
% AnovaData = [];
% AnovaData = [PlotData{1}(AnovaPCIdx,AnovaIdx); PlotData{2}(AnovaPCIdx,AnovaIdx); PlotData{3}(AnovaPCIdx,AnovaIdx); PlotData{4}(AnovaPCIdx,AnovaIdx)]';
% 
% [p,origtbl,stats]  = anova1(AnovaData)
% 
% [results,~,~,gnames] = multcompare(stats,'CriticalValueType','lsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
%%
clearvars Avg Sem Order;
Order = [3 4 1 2];
Case = 0;
for i = 1:2
    for j = 1:2
        Case = Case + 1;
        Avg{Order(Case),1} = flip(nanmean(PD{i,j},2),1);
        Sem{Order(Case),1} = flip(nanstd(PD{i,j},0,2)./sqrt(size(PD{i,j},2)),1);
    end
end

Color = [197 197 197; 138 138 138; 0 0 0]./255;

fig = figure('Position',[50 50 125 110]);
hold on
for j = 1:4
    for i = 1:3
        bar(j-0.3+(i-1)*0.3,Avg{j}(i),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.3+(i-1)*0.3,Avg{j}(i),Sem{j}(i),Sem{j}(i),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([0.4 4.6]);
xticks([1:4]);
xticklabels({'AG','ANG','VG','VNG'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Projection distance (0 to 1s)','FontName','Arial','FontSize',6);
ylim([0 15]);
    
mkdir('Figure'); cd('Figure');
saveas(fig,['Fig 2h, Unisensory Projection Distance .svg']);
cd ../
%% ANOVA test 
% clearvars Color Avg Sem
% 
% clearvars StimType LocoType
% 
% StimType{1,1} = 'Vg'; StimType{1,2} = 'Vng'; StimType{2,1} = 'Ag'; StimType{2,2} = 'Ang';
% LocoType{1} = 'KM'; LocoType{2} = 'SM'; LocoType{3} = 'S';
% 
% Data = [];
% Type1 = [];
% Type2 = [];
% N = 0;
% for i = 1:2
%     for j = 1:2
%         for k = 1:3
%             for n = 1:size(PD{i,j},2)
%                 N = N + 1;
%                 Data(N,1) = PD{i,j}(k,n);
%                 Type1{N,1} = StimType{i,j};
%                 Type2{N,1} = LocoType{k};
%             end
%         end
%     end
% end
% 
% [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
%     "Varnames",["Stim","Loco"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','lsd');
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
% 
% clearvars StimType LocoType Data Type1 Type2 N

%%

%CalBin = DistBin;
CalBin = (Pre+DistTime2)*ImgHz;

clearvars Avg Sem Order;
Order = [2 1];
for i = 1:2
    for j = 1:3
        GNG_D_Avg{i}(j,:) = nanmean(GNG_D{i,1}{j}(CalBin,:),1);
        Avg(Order(i),j) = nanmean(nanmean(GNG_D{i,1}{j}(CalBin,:),1),2);
        Sem(Order(i),j) = nanstd(nanmean(GNG_D{i,1}{j}(CalBin,:),1),0,2)./sqrt(size(GNG_D{i,1},2));
    end
end

Avg = flip(Avg,2); Sem = flip(Sem,2);

Color = [197 197 197; 138 138 138; 0 0 0]./255;

fig = figure('Position',[50 50 100 110]);
hold on
for i = 1:3
    for j = 1:2
        bar(j-0.275+(i-1)*0.275,Avg(j,i),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.275+(i-1)*0.275,Avg(j,i),Sem(j,i),Sem(j,i),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([0.4 2.6]);
xticks([1:2]);
xticklabels({'Aud','Vis'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Euclidean distance','FontName','Arial','FontSize',6);
ylim([0 25]);

    
mkdir('Figure'); cd('Figure');
saveas(fig,['Fig 2i, Unisensory GNG Distance .svg']);
cd ../


%clearvars StimType LocoType Data Type1 Type2 N

%% ANOVA test
% StimType{1,1} = 'Vis'; StimType{2,1} = 'Aud';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for i = 1:2
%     for k = 1:3
%         for n = 1:size(GNG_D{i,1}{k},2)
%             N = N + 1;
%             Data(N,1) = nanmean(GNG_D{i,1}{k}(CalBin,n),1);
%             Type1{N,1} = StimType{i};
%             Type2{N,1} = LocoType{k};
%         end
%     end
% end
% 
% [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
%     "Varnames",["Stim","Loco"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','lsd');
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))


%%

%CalBin = DistBin;
CalBin = (Pre+DistTime2)*ImgHz;

clearvars Avg Sem;

for j = 1:2
    ConPDAvg(:,j) = nanmean(PD{3,j},2); ConPDSem(:,j) = nanstd(PD{3,j},0,2)./sqrt(size(PD{3,j},2));
end

for i = 1
    for j = 1:2
        Avg{j,1} = flip(nanmean(PD{i,j},2),1);
        Sem{j,1} = flip(nanstd(PD{i,j},0,2)./sqrt(size(PD{i,j},2)),1);
    end
end

for j = 1:4
    InconPDAvg(:,j) = nanmean(PD{4,j},2); InconPDSem(:,j) = nanstd(PD{4,j},0,2)./sqrt(size({4,j},2));
end

ConPDAvg = flip(ConPDAvg,1); ConPDSem = flip(ConPDSem,1);
InconPDAvg = flip(InconPDAvg,1); InconPDSem = flip(InconPDSem,1);
%
Color = [197 197 197; 138 138 138; 0 0 0]./255;

fig = figure('Position',[50 50 125 110]);
hold on

for j = 1:2
    for i = 1:3
        bar(j-0.3+(i-1)*0.3,ConPDAvg(i,j),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.3+(i-1)*0.3,ConPDAvg(i,j),ConPDSem(i,j),ConPDSem(i,j),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

for j = 1:2
    for i = 1:3
        bar(2+j-0.3+(i-1)*0.3,Avg{j}(i),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(2+j-0.3+(i-1)*0.3,Avg{j}(i),Sem{j}(i),Sem{j}(i),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([0.4 4.6]);
xticks([1:4]);
xticklabels({'ConG','ConNG','VG','VNG'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Projection distance (0 to 1s)','FontName','Arial','FontSize',6);
ylim([0 15]);

mkdir('Figure'); cd('Figure');
saveas(fig,['Fig S3a, Multisensory Con Projection Distance .svg']);
cd ../
% 
% AnovaData = [];
% 
% AnovaData = [PD{1,1}' PD{1,2}' PD{3,1}' PD{3,2}'];
% anova1(AnovaData)

fig = figure('Position',[50 50 125 110]);
hold on
Order = [4 3 2 1];
for j = 1:4
    for i = 1:3
        bar(j-0.275+(i-1)*0.275,InconPDAvg(i,Order(j)),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.275+(i-1)*0.275,InconPDAvg(i,Order(j)),InconPDSem(i,Order(j)),InconPDSem(i,Order(j)),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([0.4 4.6]);
xticks([1:4]);
xticklabels({'V dom','A dom','V dom','A dom'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Projection distance (0 to 1s)','FontName','Arial','FontSize',6);
ylim([0 20]);

mkdir('Figure'); cd('Figure');
saveas(fig,['Fig S3e, Multisensory Incon Projection Distance.svg']);
cd ../

%% ANOVA test
% clearvars StimType LocoType Data Type1 Type2 N
% 
% StimType{1,1} = '[A dom] Ig'; StimType{2,1} = '[V dom] Ig'; StimType{3,1} = '[A dom] Ing'; StimType{4,1} = '[V dom] Ing';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for i = 1:4
%     for k = 1:3
%         for n = 1:size(PD{4,i},2)
%             N = N + 1;
%             Data(N,1) = PD{4,i}(k,n);
%             Type1{N,1} = StimType{i};
%             Type2{N,1} = LocoType{k};
%         end
%     end
% end
% 
% [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
%     "Varnames",["Type","Loco"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','lsd');
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
% 
% for i = 1:66
%     if strcmp(tbl{i,1},'Type=[A dom] Ig,Loco=S') == 1
%         if strcmp(tbl{i,2},'Type=[V dom] Ig,Loco=S') == 1
%             tbl(i,3)
%         end
%     elseif strcmp(tbl{i,1},'Type=[V dom] Ig,Loco=S') == 1
%         if strcmp(tbl{i,2},'Type=[A dom] Ig,Loco=S') == 1
%             tbl(i,3)
%         end
%     end
% end
% 
% clearvars Color
%%
close all;

%CalBin = DistBin;
CalBin = (Pre+DistTime2)*ImgHz;

Color = [197 197 197; 138 138 138; 0 0 0]./255;

for j = 1:3
    GNGAvg(j,1) = nanmean(nanmean(GNG_D{3,1}{j}(CalBin,:),1)); GNGSem(j,1) = nanstd(nanmean(GNG_D{3,1}{j}(CalBin,:),1))./sqrt(nS);
    GNGAvg(j,2) = nanmean(nanmean(GNG_D{4,1}{j}(CalBin,:),1)); GNGSem(j,2) = nanstd(nanmean(GNG_D{4,1}{j}(CalBin,:),1))./sqrt(nS);
    GNGAvg(j,3) = nanmean(nanmean(GNG_D{4,2}{j}(CalBin,:),1)); GNGSem(j,3) = nanstd(nanmean(GNG_D{4,2}{j}(CalBin,:),1))./sqrt(nS);
end

GNGAvg = flip(GNGAvg,1);
GNGSem = flip(GNGSem,1);

clearvars Avg Sem Order;
Order = [1 1];
for i = 1
    for j = 1:3
        GNG_D_Avg{i}(j,:) = nanmean(GNG_D{i,1}{j}(CalBin,:),1);
        Avg(Order(i),j) = nanmean(nanmean(GNG_D{i,1}{j}(CalBin,:),1),2);
        Sem(Order(i),j) = nanstd(nanmean(GNG_D{i,1}{j}(CalBin,:),1),0,2)./sqrt(size(GNG_D{i,1},2));
    end
end
Avg = flip(Avg,2); Sem = flip(Sem,2);

fig = figure('Position',[50 50 100 110]);
hold on

Order = [1 3 2];
for j = 1
    for i = 1:3
        bar(j-0.275+(i-1)*0.275,GNGAvg(i,Order(j)),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.275+(i-1)*0.275,GNGAvg(i,Order(j)),GNGSem(i,Order(j)),GNGSem(i,Order(j)),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

for i = 1:3
    for j = 1
        bar(1+j-0.275+(i-1)*0.275,Avg(j,i),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(1+j-0.275+(i-1)*0.275,Avg(j,i),Sem(j,i),Sem(j,i),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([0.4 2.6]);
xticks([1:2]);
xticklabels({'Con','Vis'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Euclidean distance at 1s','FontName','Arial','FontSize',6);
ylim([0 25]);

mkdir('Figure'); cd('Figure');
saveas(fig,['Fig S3b, Multisensory Con GNG Distance .svg']);
cd ../
% 
% AnovaData = [];
% 
% AnovaData = [GNG_D{1,1}{1}(CalBin,:)' GNG_D{1,1}{2}(CalBin,:)' GNG_D{1,1}{3}(CalBin,:)' GNG_D{3,1}{1}(CalBin,:)' GNG_D{3,1}{2}(CalBin,:)' GNG_D{3,1}{3}(CalBin,:)'];
% anova1(AnovaData)
% 
%
fig = figure('Position',[50 50 100 110]);
hold on

Order = [1 3 2];
for j = 2:3
    for i = 1:3
        bar(j-0.275+(i-1)*0.275,GNGAvg(i,Order(j)),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.275+(i-1)*0.275,GNGAvg(i,Order(j)),GNGSem(i,Order(j)),GNGSem(i,Order(j)),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

xlim([1.4 3.6]);
xticks([2:3]);
xticklabels({'V dom','A dom'})
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Euclidean distance at 1s','FontName','Arial','FontSize',6);
ylim([0 25]);

mkdir('Figure'); cd('Figure');
saveas(fig,['Fig S3f, Multisensory Incon GNG Distance .svg']);
cd ../

%% ANOVA test
%clearvars StimType LocoType Data Type1 Type2 N

% CalBin = (Pre+DistTime2)*ImgHz;
% 
% StimType{1,1} = ' A dom'; StimType{2,1} = 'V dom';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for i = 1:2
%     for k = 1:3
%         for n = 1:size(GNG_D{4,i}{k},2)
%             N = N + 1;
%             Data(N,1) = nanmean(GNG_D{4,i}{k}(CalBin,n),1);
%             Type1{N,1} = StimType{i};
%             Type2{N,1} = LocoType{k};
%         end
%     end
% end
% 
% [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
%     "Varnames",["Dom","Loco"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','lsd');
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))

%%
Color = [1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 1; 1 0 1; 0 0 0; 0.6 0.6 0.6];

CalBin = (Pre+DistTime2)*ImgHz;

Order = [6 1 2 3 4 5 7 8 9];
AvgUniD = nanmean(Uni_D(CalBin,:,:),1);
for State = 1:3
    if State == 1; Session = SIdx; elseif State == 2; Session = MIdx; elseif State == 3; Session = [MIdx; SIdx;]; end;
    
    Fig = figure('Position',[50 50 150 120]);
    hold on
    for i = 1:numel(Order)
        Idx = Order(i);
        
        Avg = nanmean(AvgUniD(:,Idx,Session),3);
        Sem = nanstd(AvgUniD(:,Idx,Session),0,3)/sqrt(numel(Session));
        bar(i,Avg,0.6,'FaceColor','none','EdgeColor',Color(Idx,:),'lineWidth',1,'lineStyle','-');
        errorbar(i,Avg,Sem,Sem,'CapSize',3,'lineWidth',0.75,'color',Color(Idx,:));
        
        %     X = []; X(1:numel(S_Session)) =  i;
        %     scatter(X,squeeze(AvgUniD(:,Idx,S_Session)),'MarkerEdgeColor',Color(Order(i),:),'MarkerFaceColor','none');
        %     X = []; X(1:numel(M_Session)) =  i;
        %     scatter(X,squeeze(AvgUniD(:,Idx,M_Session)),'MarkerEdgeColor',Color(Order(i),:),'MarkerFaceColor',[0.5 0.5 0.5]);
    end
    
    xticks([1:9]); xlim([0.4 9.6]);
    xticklabels({'A_G - A_N_G','V_G - V_N_G','V_G - A_G','V_G - A_N_G','V_N_G - A_G','V_N_G - A_N_G','CON','V dom','A dom'});
    xtickangle(35);
    ylim([0 25]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    ylabel('Euclidean distance at 1s','FontName','Arial','FontSize',6);
    
    mkdir('Figure'); cd('Figure');
    if State == 1;
        saveas(Fig,'Fig S2f, Stationary Unisensory Distance.svg');
    elseif State == 2
        saveas(Fig,'Fig S2f, Moving Unisensory Distance.svg');
    elseif State == 3
        saveas(Fig,'Fig S2f, Total Unisensory Distance.svg');
    end
    cd ../
end
close all;
%% ANOVA test
% AnovaData = [];
% for i = 1:6
%     AnovaData = cat(2,AnovaData,squeeze(Uni_D(CalBin,Order(i),11:20)));
% end
% 
% [p,origtbl,stats]  = anova1(AnovaData)
% 
% [results,~,~,gnames] = multcompare(stats,'CriticalValueType','lsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
%%

Color = [255 0 0; 243 152 0; 0 0 255; 0 152 243]./255;
Style = {'-','-','-','-'};

%CalBin = DistBin;
CalBin = (Pre+DistTime2)*ImgHz;;

clearvars Avg Sem Data;
for i = 1:2
    for j = 1:3
        Avg(i,j) = nanmean(nanmean(MvmtDv{3,i}{j}(CalBin,:),1),2);
        Sem(i,j) = nanstd(nanmean(MvmtDv{3,i}{j}(CalBin,:),1),0,2)./sqrt(size(MvmtDv{3,i}{j},2));
        Data{i}(j,:) = nanmean(MvmtDv{3,i}{j}(CalBin,:),1);
    end
end

for i = 1:2
    for j = 1:3
        Avg(i+2,j) = nanmean(nanmean(MvmtDa{3,i}{j}(CalBin,:),1),2);
        Sem(i+2,j) = nanstd(nanmean(MvmtDa{3,i}{j}(CalBin,:),1),0,2)./sqrt(size(MvmtDa{3,i}{j},2));
        Data{i+2}(j,:) = nanmean(MvmtDa{3,i}{j}(CalBin,:),1);
    end
end

Avg = flip(Avg,2);
Sem = flip(Sem,2);

Fig = figure('Position',[50 50 180 110]);
hold on

xAxis = [1 2 3; 1.3 2.3 3.3; 4 5 6; 4.3 5.3 6.3];

for i = 1:4
    plot(xAxis(i,:),Avg(i,:),'color',Color(i,:),'lineStyle',Style{i},'lineWidth',1);
    errorbar(xAxis(i,:),Avg(i,:),Sem(i,:),Sem(i,:),'color',Color(i,:),'lineStyle',Style{i},'lineWidth',0.75,'CapSize',3);
end

xlim([0.4 6.9]);
ylim([0 20]);

xticks([1.15 2.15 3.15 4.15 5.15 6.15]);
xticklabels({'S','SM','KM','S','SM','KM'});

set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Locomotion state','FontName','Arial','FontSize',6);
ylabel('Euclidean distance','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(Fig,'Fig 3b, Con Distance to unisensory.svg');
cd ../

%% ANOVA test
% clearvars StimType LocoType Data Type1 Type2 Type3 Type4 N
% 
% CalBin = (Pre+DistTime2)*ImgHz;
% 
% StimType{1,1} = 'Cg to Vg'; StimType{2,1} = 'Cng to Vng'; StimType{3,1} = 'Cg to Ag'; StimType{4,1} = 'Cng to Ang';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% clearvars Avg Sem Data;
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for i = 1:2
%     for j = 1:3
%         for jj = 1:size(MvmtDv{3,i}{j},2)
%             N = N+1;
%             Data(N,1) = MvmtDv{3,i}{j}(CalBin,jj);
%             Type1{N,1} = StimType{i};
%             Type2{N,1} = LocoType{j};
%         end
%     end
% end
% 
% for i = 1:2
%     for j = 1:3
%         for jj = 1:size(MvmtDa{3,i}{j},2)
%             N = N+1;
%             Data(N,1) = MvmtDa{3,i}{j}(CalBin,jj);
%             Type1{N,1} = StimType{i+2};
%             Type2{N,1} = LocoType{j};
%         end
%     end
% end
% 
% 
% [p,origtbl,stats] = anovan(Data,{Type1 Type2},"Model","interaction", ...
%     "Varnames",["Stim","Loco"]);
% 
% [results,~,~,gnames] = multcompare(stats,"Dimension",[1 2],'CriticalValueType','hsd');
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
% 

%%
%CalBin = DistBin;
CalBin = (Pre+DistTime2)*ImgHz;

Color = [243 152 0; 243 152 0; 255 0 0; 255 0 0; 0 0 255; 0 0 255; 0 152 243; 0 152 243]./255;
Style = {'-','-','-','-','-','-','-','-'};

clearvars Avg Sem;
for i = 1:4
    for j = 1:3
        Avg(i,j) = nanmean(nanmean(MvmtInconDv{i,j}(CalBin,:),1),2);
        Sem(i,j) = nanstd(nanmean(MvmtInconDv{i,j}(CalBin,:),1),0,2)./sqrt(size(MvmtInconDv{i,j},2));
    end
end

for i = 1:4
    for j = 1:3
        Avg(i+4,j) = nanmean(nanmean(MvmtInconDa{i,j}(CalBin,:),1),2);
        Sem(i+4,j) = nanstd(nanmean(MvmtInconDa{i,j}(CalBin,:),1),0,2)./sqrt(size(MvmtInconDa{i,j},2));
    end
end

Avg = flip(Avg,2);
Sem = flip(Sem,2);

Fig = figure('Position',[50 50 180 110]);
hold on

xAxis = [4.8 5.8 6.8; 1.3 2.3 3.3; 4.5 5.5 6.5; 1 2 3];
xAxis = [xAxis; xAxis+7.5];

for i = 1:8
    plot(xAxis(i,:),Avg(i,:),'color',Color(i,:),'lineStyle',Style{i},'lineWidth',1);
    errorbar(xAxis(i,:),Avg(i,:),Sem(i,:),Sem(i,:),'color',Color(i,:),'lineStyle',Style{i},'lineWidth',0.75,'CapSize',3);
end

xlim([0.4 14.8]);
ylim([0 25]);

xticks([1.15 2.15 3.15 4.65 5.65 6.65 8.65 9.65 10.65 12.15 13.15 14.15]);
xticklabels({'S','SM','KM','S','SM','KM','S','SM','KM','S','SM','KM'});

set(gca,'TickDir','out','FontName','Arial','FontSize',6);

xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Euclidean distance to V','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(Fig,'Fig 3d, Incon Mean Distance to unisensory.svg');
cd ../

%%
close all;
clearvars X Y;

CalBin = (Pre+DistTime2)*ImgHz;

Color = [0 0 0; 100 100 100; 200 200 200]./255;  
AvgColor = [255 0 0; 255 100 100; 255 180 180]./255;

for j = 1:2
    Fig = figure('Position',[50 50 90 110]);
    hold on
    
    for N = 1:nS
        for jj = 1:3
            X(N,jj) = nanmean(MvmtDv{4,j}{jj}(CalBin,N),1);
        end
        Y(N,:) = ADomRate(j,:,N).*100;
    end

    for N = 1:nS
        plot(X(N,:),Y(N,:),'lineWidth',0.25,'color',[200 200 200]./255);
        alpha(0.2);
    end
    
    for jj = 1:3
        scatter(X(:,jj),Y(:,jj),15,'MarkerEdgeColor','none','MarkerFaceColor',Color(jj,:));
        alpha(0.5);
    end
    
    AvgX = nanmean(X,1);
    AvgY = nanmean(Y,1);
    SemX = nanstd(X,0,1);%./sqrt(size(X,1));
    SemY = nanstd(Y,0,1);%./sqrt(size(Y,1));
    
     for jj = 1:3
         errorbar(AvgX(jj),AvgY(jj),SemY(jj),SemY(jj),SemX(jj),SemX(jj),'lineWidth',1.5,'color',AvgColor(jj,:),'CapSize',0);
     end
    
    xlim([0 11]);
    %xticks([0:5:15]);
    ylim([0 100]);
    set(gca,'TickDir','out','FontName','Arial','FontSize',6);
    xlabel('Euclidean distance to V','FontName','Arial','FontSize',6);
    ylabel('Aud dominance rate (%)','FontName','Arial','FontSize',6);
    
    mkdir('Figure'); cd('Figure');
    if j == 1; saveas(Fig,'Fig 3g, AV Distance A Dom rate.svg');
    elseif j == 2; saveas(Fig,'Fig 3g, VA Distance A Dom rate.svg');
    end
    cd ../
end