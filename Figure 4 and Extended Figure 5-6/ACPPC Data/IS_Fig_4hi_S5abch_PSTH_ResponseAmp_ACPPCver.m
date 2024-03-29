% This code calculates response amplitudes and makes figures of PSTH and
% response amplitudes
% Run this code in the folder where the exmaple data are.

clear variables; close all; clc;

MatFiles = FindMatFiles();

nS = 0; nM = 0;
for MatIdx = 1:numel(MatFiles)
    if strcmp(MatFiles{MatIdx}(1:6),'Moving') == 1
        nM = nM + 1;
        M_Mat{nM,1} = MatFiles{MatIdx};
    elseif strcmp(MatFiles{MatIdx}(1:10),'Stationary') == 1
        nS = nS + 1;
        S_Mat{nS,1} = MatFiles{MatIdx};
    end
end

for j = 1:4
    InconF{1,j} = []; %Stationary
    InconF{2,j} = []; %Moving
    for i = 1:3
        InconAmp{i,j} = [];
    end
end

for i = 1:3
    for ii = 1:3
        StimF{i}{ii,1} = [];
        StimF{i}{ii,3} = [];
        StimAmp{i}{ii,1} = [];
        StimAmp{i}{ii,3} = [];
    end
end

for MatIdx = 1:nS
    load(S_Mat{MatIdx});

    Post = 1;
    Pre = 0.2;

    nPre = PSTHPre*ImgHz;

    for i = 1:3
        for StimIdx = [1 3]
            StimF{3}{i,StimIdx} = cat(2,StimF{3}{i,StimIdx},nanmean(EventF{i,StimIdx},3));
            StimAmp{3}{i,StimIdx} = cat(2,StimAmp{3}{i,StimIdx},nanmean(nanmean(EventF{i,StimIdx}(nPre+1:(PSTHPre+Post)*ImgHz,:,:),1) - nanmean(EventF{i,StimIdx}((PSTHPre-Pre)*ImgHz+1:nPre,:,:),1),3));
        end
    end

    for j = 1:4
        InconF{1,j} = cat(2,InconF{1,j},nanmean(EventF{4,j},3));
        InconAmp{3,j} = cat(2,InconAmp{3,j},nanmean(nanmean(EventF{4,j}(nPre+1:(PSTHPre+Post)*ImgHz,:,:),1) - nanmean(EventF{4,j}((PSTHPre-Pre)*ImgHz+1:nPre,:,:),1),3));
    end
end

for MatIdx = 1:nM
    load(M_Mat{MatIdx});

    for i = 1:3
        for StimIdx = [1 3]
            KMIdx = []; SMIdx = [];
            for k = 1:size(EventSpeed{i,StimIdx},1)
                if mean(EventSpeed{i,StimIdx}(k,1:KMTime*EventHz)) > KMThreshold
                    KMIdx = cat(1,KMIdx,k);
                else
                    SMIdx = cat(1,SMIdx,k);
                end
            end
            StimF{2}{i,StimIdx} = cat(2,StimF{2}{i,StimIdx},nanmean(EventF{i,StimIdx}(:,:,SMIdx),3));
            StimF{1}{i,StimIdx} = cat(2,StimF{1}{i,StimIdx},nanmean(EventF{i,StimIdx}(:,:,KMIdx),3));
            StimAmp{2}{i,StimIdx} = cat(2,StimAmp{2}{i,StimIdx},nanmean(nanmean(EventF{i,StimIdx}(nPre+1:(PSTHPre+Post)*ImgHz,:,SMIdx),1) - nanmean(EventF{i,StimIdx}((PSTHPre-Pre)*ImgHz+1:nPre,:,SMIdx),1),3));
            StimAmp{1}{i,StimIdx} = cat(2,StimAmp{1}{i,StimIdx},nanmean(nanmean(EventF{i,StimIdx}(nPre+1:(PSTHPre+Post)*ImgHz,:,KMIdx),1) - nanmean(EventF{i,StimIdx}((PSTHPre-Pre)*ImgHz+1:nPre,:,KMIdx),1),3));
        end
    end

    for j = 1:4
        InconF{2,j} = cat(2,InconF{2,j},nanmean(EventF{4,j},3));
        KMIdx = []; SMIdx = [];
           for k = 1:size(EventSpeed{4,j},1)
                if mean(EventSpeed{4,j}(k,1:KMTime*EventHz)) > KMThreshold
                    KMIdx = cat(1,KMIdx,k);
                else
                    SMIdx = cat(1,SMIdx,k);
                end
           end
           InconAmp{1,j} = cat(2,InconAmp{1,j},nanmean(nanmean(EventF{4,j}(nPre+1:(PSTHPre+Post)*ImgHz,:,KMIdx),1) - nanmean(EventF{4,j}((PSTHPre-Pre)*ImgHz+1:nPre,:,KMIdx),1),3));
           InconAmp{2,j} = cat(2,InconAmp{2,j},nanmean(nanmean(EventF{4,j}(nPre+1:(PSTHPre+Post)*ImgHz,:,SMIdx),1) - nanmean(EventF{4,j}((PSTHPre-Pre)*ImgHz+1:nPre,:,SMIdx),1),3));
    end
end

%%

xAxis = -PSTHPre+1/ImgHz:1/ImgHz:PSTHPost;

Color = [0 0 0; 120 120 120; 180 180 180]./255;

YLim = [-0.5 1.5];

Style{1} = '-'; Style{3} = '-';

StateName{1} = 'Stationary'; StateName{2} = 'Moving';
StimName{1,1} = 'Fig S5b, Vgo'; StimName{1,3} = 'Fig S5b, Vnogo'; StimName{2,1} = 'Fig S5a, Ago'; StimName{2,3} = 'Fig S5a, Anogo';

for i = 1:2
    for j = [1 3]
        fig = figure('Position',[50 50 60 110]);
        hold on

        for k = 1:3
            Avg1 = nanmean(StimF{k}{i,j},2);
            Sem1 = nanstd(StimF{k}{i,j},0,2)./sqrt(size(StimF{k}{i,j},2));

            errorshade(xAxis,Avg1',Avg1'+Sem1',Avg1'-Sem1',Color(k,:));
            alpha(0.4);
        end
        for k = 1:3
            Avg1 = nanmean(StimF{k}{i,j},2);
            Sem1 = nanstd(StimF{k}{i,j},0,2)./sqrt(size(StimF{k}{i,j},2));

            plot(xAxis,Avg1,'lineWidth',1,'color',Color(k,:));
        end
        plot(xAxis,Avg1,'lineWidth',1,'color',Color(k,:));

        set(gca,'TickDir','out','FontName','Arial','FontSize',6);
        xlabel('Time (s)','FontName','Arial','FontSize',6);
        ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

        line([0 0],YLim,'lineWidth',0.5,'color',[160 160 160]./255);
        line([1 1],YLim,'lineWidth',0.5,'color',[160 160 160]./255,'lineStyle',':');

        xlim([-1 2]); yticks([0 1 2]);
        ylim([YLim]);

        mkdir('Figure'); cd('Figure');
        saveas(fig,[StimName{i,j} ' PSTH.svg']);
        cd ../
    end
end


Color = [0 0 0; 138 138 138; 197 197 197;]./255;

Center = [3 4 1 2];
Idx = 0;
fig = figure('Position',[50 50 145 110]);
hold on
for i = 1:2
    for j = [1 3]
        Idx = Idx + 1;
        for k = 1:3
            Avg = nanmean(StimAmp{k}{i,j},2);
            Sem = nanstd(StimAmp{k}{i,j},0,2)./sqrt(size(StimAmp{k}{i,j},2));

            bar(Center(Idx)+0.3-(k-1)*0.3,Avg,0.25,'lineWidth',1,'EdgeColor',Color(k,:),'FaceColor','none');
            errorbar(Center(Idx)+0.3-(k-1)*0.3,Avg,Sem,Sem,'CapSize',3,'color',Color(k,:),'lineWidth',0.75)

        end
    end
end
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
xlabel('Type','FontName','Arial','FontSize',6);
ylabel('Response amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);

xticks([1:4]); xticklabels({'AG','ANG','VG','VNG'});
xlim([0.4 4.6]); %yticks([0 1 2]);
ylim([0 0.6]);
%ylim([]);
mkdir('Figure'); cd('Figure');
saveas(fig,['Fig S5c, Unisensory response amplitude.svg']);
cd ../

clearvars Avg Sem Data;
% 
% clearvars StimType LocoType Data Type1 Type2 Type3 Type4 N
% 
% StimType{1,1} = 'Vg'; StimType{1,3} = 'Vng'; StimType{2,1} = 'Ag'; StimType{2,3} = 'Ang';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% 
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for k = 1:3
%     for i = [1 2]
%         for j = [1 3]
%             for jj = 1:size(StimAmp{k}{i,j},2)
%                 N = N+1;
%                 Data(N,1) = StimAmp{k}{i,j}(jj);
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
% %[results,~,~,gnames] = multcompare(stats,"Dimension",[2],'CriticalValueType','hsd');
% 
% tbl = array2table(results(:,[1 2 6]),"VariableNames", ...
%     ["Group A","Group B","P-value"]);
% tbl.("Group A")=gnames(tbl.("Group A"));
% tbl.("Group B")=gnames(tbl.("Group B"))
%%
clearvars StateName StimName

xAxis = -PSTHPre+1/ImgHz:1/ImgHz:PSTHPost;

Color = [0.6 0.6 0.6; 0 0 0];

Style{1} = '-'; Style{3} = '-';

StateName{1} = 'Fig 4h, Stationary'; StateName{2} = 'Fig 4i, Moving';
StimName{1} = 'AgoVnogo'; StimName{3} = 'AnogoVgo';

for i = 1:2
    for j = [1 3]
        for ii = 1:size(InconF{i,j},1)
            p(ii,1) =  signrank(InconF{i,j}(ii,:),InconF{i,j+1}(ii,:));
        end

        MildIdx = find(p < 0.05 & p >= 0.01);
        StrongIdx = find(p < 0.01);

        fig = figure('Position',[50 50 60 110]);
        hold on
        Avg1 = nanmean(InconF{i,j},2);
        Sem1 = nanstd(InconF{i,j},0,2)./sqrt(size(InconF{i,j},2));

        Avg2 = nanmean(InconF{i,j+1},2);
        Sem2 = nanstd(InconF{i,j+1},0,2)./sqrt(size(InconF{i,j+1},2));

        errorshade(xAxis,Avg1',Avg1'+Sem1',Avg1'-Sem1',Color(1,:));
        alpha(0.4);
        errorshade(xAxis,Avg2',Avg2'+Sem2',Avg2'-Sem2',Color(2,:));
        alpha(0.4);

        plot(xAxis,Avg1,'lineWidth',1,'color',Color(1,:),'lineStyle',Style{j});
        plot(xAxis,Avg2,'lineWidth',1,'color',Color(2,:),'lineStyle',Style{j});

        % for ii = MildIdx';
        %     line([xAxis(ii)-(1/ImgHz)/2 xAxis(ii)+(1/ImgHz)/2],[YLim(2) YLim(2)],'color',[0.6 0.6 0.6],'lineWidth',3);
        % end
        % for ii = StrongIdx';
        %     line([xAxis(ii)-(1/ImgHz)/2 xAxis(ii)+(1/ImgHz)/2],[YLim(2) YLim(2)],'color',[0 0 0],'lineWidth',3);
        % end

        set(gca,'TickDir','out','FontName','Arial','FontSize',6);
        xlabel('Time (s)','FontName','Arial','FontSize',6);
        ylabel('ΔF/F (%)','FontName','Arial','FontSize',6);

        line([0 0],YLim,'lineWidth',0.5,'color',[160 160 160]./255);
        line([1 1],YLim,'lineWidth',0.5,'color',[160 160 160]./255,'lineStyle',':');

        xlim([-1 2]); yticks([0 1 2]);
        ylim([YLim]);

        mkdir('Figure'); cd('Figure');
        saveas(fig,[StateName{i} ' ' StimName{j} ' PSTH.svg']);
        cd ../
    end
end

clearvars Avg Sem Style;

fig = figure('Position',[50 50 100 110]);
hold on
xTicks = [5.1 6.1 7.1; 1.7 2.7 3.7; 4.7 5.7 6.7; 1.3 2.3 3.3];
Style = {':',':','-','-'};

for i = 1:3
    for j = 1:4
        Avg(i,j) = nanmean(InconAmp{i,j});
        Sem(i,j) = nanstd(InconAmp{i,j})./sqrt(numel(InconAmp{i,j}));
    end
end

Avg = flip(Avg,1);
Sem = flip(Sem,1);
for j = 1:4
    errorbar(xTicks(j,:),Avg(:,j),Sem(:,j),Sem(:,j),'CapSize',3,'color','k','lineWidth',.75,'lineStyle',Style{j});
end


xlim([0.7 7.7]); ylim([0 0.7]);
xticks([1.5 2.5 3.5 4.9 5.9 6.9]);
xticklabels({'S','SM','KM','S','SM','KM'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6);
ylabel('Repsonse amplitude, ΔF/F (%)','FontName','Arial','FontSize',6);
xlabel('V dom          A dom','FontName','Arial','FontSize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig S5h, Incon response amplitude.svg');
cd ../
%%
% clearvars StimType LocoType Data Type1 Type2 Type3 Type4 N
% 
% StimType{1,1} = '[A dom] Ig'; StimType{2,1} = '[V dom] Ig'; StimType{3,1} = '[A dom] Ing'; StimType{4,1} = '[V dom] Ing';
% LocoType{1,1} = 'KM'; LocoType{2,1} = 'SM'; LocoType{3,1} = 'S';
% 
% clearvars Avg Sem Data;
% Data = []; Type1 = []; Type2 = [];
% N = 0;
% for i = 1:3
%     for j = 1:4
%         for jj = 1:size(InconAmp{i,j},2)
%             N = N+1;
%             Data(N,1) = InconAmp{i,j}(jj);
%             Type1{N,1} = StimType{j};
%             Type2{N,1} = LocoType{i};
%         end
%     end
% end
% 
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