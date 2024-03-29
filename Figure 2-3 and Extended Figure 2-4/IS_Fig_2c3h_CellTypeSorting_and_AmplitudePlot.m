% This code sort neuronal response type based on activity change by unsensory stimuli
% Cell type sorting : Significant based; For both G/NG responsive neuron,
% classified into G or NG according to response amplitude.
% In other words, V and A neurons have statistically stimilar G/NG responses.
% Run this code where the example data are.

% Variable named "nCell" shows number of each response type in 
% stationary (1st column) and moving (2nd column) sessions.
% 1st row: Vgo preferring neuron; 
% 2nd row: Vnogo preferring neuron; 
% 3rd row: both Vgo and Vnogo preferring neuron (Vgo&nogo neuron);
% 4th row: Ago preferring neuron; 
% 5th row: Anogo preferring neuron; 
% 6th row: both Ago and Anogo preferring neuron (Ago&nogo neuron);
% 7th row: Dual responsive (one visual; one auditory)
% 8th row: Triple responsive (two visual & one auditory and vice versa)
% 9th row: all responsive neuron
% 10th row: non-responsive to all the unisensory stimuli

clear variables; close all; clc;
%set(0,'DefaultFigureVisible','on')

Mat = FindMatFiles();

Target = [1 1; 1 3; 2 1; 2 3; 3 1; 3 3; 4 1; 4 2; 4 3; 4 4;];
Pre = 0.2; %Time epoch for pre-stimulus activity (in second)
Post = 1; %Time epoch for post-stimulus activity (in second)
PValue = 0.01; %P value threshold for significant response
KMTime = 1.5;
KMThreshold = 3;

OnlyInc = 1 % 0 = calculate both increasing and decreasing neurons for Figure 3h; 1 = calculate increasing neurons only
Abs = 1 % 0 = just average increasing and decreasing activity for Figure 3h; 1 = take absolute value of respons amplitude
OnlySelective = 1 % 0 = Also consider Vgo&nogo and Ago&nogo for Figure 3h; 1 = Consider go-selective or nogo-selective neurons only

MergedAmp = [];
MergedIndAmp = [];
MergedSIndAmp = []; MergedSMIndAmp = []; MergedKMIndAmp = [];
MergedIndAmp = [];
S_H = []; M_H = []; M_IndAmp = []; S_IndAmp = []; M_AvgAmp = []; S_AvgAmp = [];
%
for i = 1:size(Target,1)
    MergedF{i,1} = [];
end
CellIdx = 0;
for MatIdx = 1:numel(Mat)
    load(Mat{MatIdx});
    %ImgHz = imgHz;
    TempAvgAmp = [];
    nCell = size(EventF{1,1},2);
    clearvars TempIndAmp IndSAmp IndSMAmp IndKMAmp

    for i = 1:size(Target,1)
        x = Target(i,1); y = Target(i,2); Amp = []; KMAmp = []; SMAmp = [];
        for k = 1:size(EventF{x,y},3)
            Amp(k,:) = nanmean(EventF{x,y}(PSTHPre*ImgHz+1:(PSTHPre+Post)*ImgHz,:,k),1) - nanmean(EventF{x,y}((PSTHPre-Pre)*ImgHz+1:PSTHPre*ImgHz,:,k),1);
            if strcmp(State,'m') == 1
                if mean(EventSpeed{x,y}(k,1:KMTime*EventHz)) > KMThreshold
                    KMAmp = cat(1,KMAmp,Amp(k,:));
                else
                    SMAmp = cat(1,SMAmp,Amp(k,:));
                end
            end
        end
        MergedF{i,1} = cat(2,MergedF{i,1},nanmean(EventF{x,y},3));
        
        for j = 1:nCell
            TempIndAmp{i,j} = Amp(:,j);
            if strcmp(State,'m') == 1
                IndSMAmp{i,j} = SMAmp(:,j);
                IndKMAmp{i,j} = KMAmp(:,j);
            elseif strcmp(State,'s') == 1
                IndSAmp{i,j} = Amp(:,j);
            end
        end
        
        AvgAmp = nanmean(Amp,1);
        TempAvgAmp = cat(1,TempAvgAmp,AvgAmp);
        for j = 1:size(Amp,2)
            p = signrank(Amp(:,j));
            if p > PValue;
                H{MatIdx}(i,j) = 0;
            elseif p < PValue & nanmean(Amp(:,j)) > 0
                H{MatIdx}(i,j) = 1;
            elseif p < PValue & nanmean(Amp(:,j)) < 0
                H{MatIdx}(i,j) = -1;
            end
        end
    end
    MergedAmp = cat(2,MergedAmp,TempAvgAmp);
    MergedIndAmp = [MergedIndAmp TempIndAmp];
    if strcmp(State,'m') == 1
        MergedSMIndAmp = [MergedSMIndAmp IndSMAmp];
        MergedKMIndAmp = [MergedKMIndAmp IndKMAmp];
        M_IndAmp = [M_IndAmp TempIndAmp];
        M_AvgAmp = [M_AvgAmp TempAvgAmp];
        M_H = cat(2,M_H,H{MatIdx});
    elseif strcmp(State,'s') == 1
        MergedSIndAmp = [MergedSIndAmp IndSAmp];
        S_IndAmp = [S_IndAmp TempIndAmp];
        S_AvgAmp = [S_AvgAmp TempAvgAmp];
        S_H = cat(2,S_H,H{MatIdx});
    end
end

%%
clearvars -except S_H M_H M_IndAmp S_IndAmp M_AvgAmp S_AvgAmp MergedSIndAmp MergedSMIndAmp MergedKMIndAmp Change OnlyInc Abs OnlySelective

H_Data = S_H; AmpData = S_IndAmp;

NonList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);

VgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
VngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
AgList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
AngList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);

VList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
AList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);
VgAgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
VgAngList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
VngAngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
VngAgList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);

V_AgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
V_AngList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
A_VgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);
A_VngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);

AllList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);

V_VgList = []; V_VngList = [];
for jj = 1:numel(VList)
    p = ranksum(AmpData{1,VList(jj)},AmpData{2,VList(jj)});
    if p < 0.05
        if abs(nanmean(AmpData{1,VList(jj)})) > abs(nanmean(AmpData{2,VList(jj)}))
            V_VgList = cat(2,V_VgList,VList(jj));
        else
            V_VngList = cat(2,V_VngList,VList(jj));
        end
    end
end

Common = intersect(VList,V_VgList);
VList = setxor(VList,Common); VgList = cat(2,VgList,V_VgList);
Common = intersect(VList,V_VngList);
VList = setxor(VList,Common); VngList = cat(2,VngList,V_VngList);
if size(VList,1) > 1
    VList = VList';
end

A_AgList = []; A_AngList = [];
for jj = 1:numel(AList)
    p = ranksum(AmpData{3,AList(jj)},AmpData{4,AList(jj)});
    if p < 0.05
        if abs(nanmean(AmpData{3,AList(jj)})) > abs(nanmean(AmpData{4,AList(jj)}))
            A_AgList = [A_AgList AList(jj)];
        else
            A_AngList = [A_AngList AList(jj)];
        end
    end
end

Common = intersect(AList,A_AgList);
AList = setxor(AList,Common); AgList = cat(2,AgList,A_AgList);
Common = intersect(AList,A_AngList);
AList = setxor(AList,Common); AngList = cat(2,AngList,A_AngList);
if size(AList,1) > 1
    AList = AList';
end

nCell(1,1) = numel(VgList);
nCell(2,1) = numel(VngList);
nCell(3,1) = numel(VList);
nCell(4,1) = numel(AgList);
nCell(5,1) = numel(AngList);
nCell(6,1) = numel(AList);
nCell(7,1) = numel(VgAgList)+ numel(VgAngList)+ numel(VngAngList)+ numel(VngAgList);
nCell(8,1) = numel(V_AgList) + numel(V_AngList)+ numel(A_VgList)+ numel(A_VngList);
nCell(9,1) = numel(AllList);
nCell(10,1) = numel(NonList);

nDual(1,1) = numel(VgAngList); nDual(2,1) = numel(VgAgList); nDual(3,1) = numel(VngAngList); nDual(4,1) = numel(VngAgList);
nTriple(1,1) = numel(V_AgList); nTriple(2,1) = numel(V_AngList); nTriple(3,1) = numel(A_VgList); nTriple(4,1) = numel(A_VngList);

%%


% BiasedCellIdx{1} = VgList;
% BiasedCellIdx{2} = VngList;
% BiasedCellIdx{3} = AgList;
% BiasedCellIdx{4} = AngList;
if OnlyInc == 1
    if OnlySelective == 1
        BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) == 1))];
        BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) == 1))];
        BiasedCellIdx{3} = [AgList(find(H_Data(3,AgList) == 1))];
        BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) == 1))];
    else
        BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) == 1)) VList(find(H_Data(1,VList) == 1))];
        BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) == 1)) VList(find(H_Data(2,VList) == 1))];
        BiasedCellIdx{3} = [AgList(find(H_Data(3,AgList) == 1)) AList(find(H_Data(3,AList) == 1))];
        BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) == 1)) AList(find(H_Data(4,AList) == 1))];
    end
elseif OnlyInc == 0
    BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) ~= 0)) VList(find(H_Data(1,VList) ~= 0))];
    BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) ~= 0)) VList(find(H_Data(2,VList) ~= 0))];
    BiasedCellIdx{3} = [AgLiaqst(find(H_Data(3,AgList) ~= 0)) AList(find(H_Data(3,AList) ~= 0))];
    BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) ~= 0))  AList(find(H_Data(4,AList) ~= 0))];
end


C_Target = [5 6 5 6];
I_Target = [9 10; 7 8; 7 8; 9 10];

for i = 1:numel(BiasedCellIdx)
    clearvars Uni UniAmp
    Idx = BiasedCellIdx{i};
    for jj = 1:numel(Idx)
        Uni = nanmean(MergedSIndAmp{i,Idx(jj)});
        Con = nanmean(MergedSIndAmp{C_Target(i),Idx(jj)});
        Incon = nanmean([MergedSIndAmp{I_Target(i,1),Idx(jj)}; MergedSIndAmp{I_Target(i,2),Idx(jj)}]);
        
         % U_Data{1,i}(jj,1) = Uni;
         % C_Data{1,i}(jj,1) = Con;
         % I_Data{1,i}(jj,1) = Incon;

         if Abs == 1
             U_Data{1,i}(jj,1) = abs(Uni);
             C_Data{1,i}(jj,1) = abs(Con);
             I_Data{1,i}(jj,1) = abs(Incon);
         elseif Abs == 0
             U_Data{1,i}(jj,1) = (Uni);
             C_Data{1,i}(jj,1) = (Con);
             I_Data{1,i}(jj,1) = (Incon);
         end
    end
end

%%

clearvars -except S_H M_H M_IndAmp S_IndAmp M_AvgAmp S_AvgAmp MergedSIndAmp MergedSMIndAmp MergedKMIndAmp U_Data C_Data I_Data nCell OnlyInc Abs OnlySelective

H_Data = M_H; AmpData = M_IndAmp;

NonList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);

VgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
VngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
AgList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
AngList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);

VList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) == 0);
AList = find(H_Data(1,:) == 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);
VgAgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
VgAngList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
VngAngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
VngAgList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);

V_AgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) == 0);
V_AngList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) == 0 & H_Data(4,:) ~= 0);
A_VgList = find(H_Data(1,:) ~= 0 & H_Data(2,:) == 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);
A_VngList = find(H_Data(1,:) == 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);

AllList = find(H_Data(1,:) ~= 0 & H_Data(2,:) ~= 0 &H_Data(3,:) ~= 0 & H_Data(4,:) ~= 0);

V_VgList = []; V_VngList = [];
for jj = 1:numel(VList)
    p = ranksum(AmpData{1,VList(jj)},AmpData{2,VList(jj)});
    if p < 0.05
        if abs(nanmean(AmpData{1,VList(jj)})) > abs(nanmean(AmpData{2,VList(jj)}))
            V_VgList = cat(2,V_VgList,VList(jj));
        else
            V_VngList = cat(2,V_VngList,VList(jj));
        end
    end
end

Common = intersect(VList,V_VgList);
VList = setxor(VList,Common); VgList = cat(2,VgList,V_VgList);
Common = intersect(VList,V_VngList);
VList = setxor(VList,Common); VngList = cat(2,VngList,V_VngList);
if size(VList,1) > 1
    VList = VList';
end

A_AgList = []; A_AngList = [];
for jj = 1:numel(AList)
    p = ranksum(AmpData{3,AList(jj)},AmpData{4,AList(jj)});
    if p < 0.05
        if abs(nanmean(AmpData{3,AList(jj)})) > abs(nanmean(AmpData{4,AList(jj)}))
            A_AgList = [A_AgList AList(jj)];
        else
            A_AngList = [A_AngList AList(jj)];
        end
    end
end

Common = intersect(AList,A_AgList);
AList = setxor(AList,Common); AgList = cat(2,AgList,A_AgList);
Common = intersect(AList,A_AngList);
AList = setxor(AList,Common); AngList = cat(2,AngList,A_AngList);
if size(AList,1) > 1
    AList = AList';
end

nCell(1,2) = numel(VgList);
nCell(2,2) = numel(VngList);
nCell(3,2) = numel(VList);
nCell(4,2) = numel(AgList);
nCell(5,2) = numel(AngList);
nCell(6,2) = numel(AList);
nCell(7,2) = numel(VgAgList)+ numel(VgAngList)+ numel(VngAngList)+ numel(VngAgList);
nCell(8,2) = numel(V_AgList) + numel(V_AngList)+ numel(A_VgList)+ numel(A_VngList);
nCell(9,2) = numel(AllList);
nCell(10,2) = numel(NonList)

nDual(1,1) = numel(VgAngList); nDual(2,1) = numel(VgAgList); nDual(3,1) = numel(VngAngList); nDual(4,1) = numel(VngAgList);
nTriple(1,1) = numel(V_AgList); nTriple(2,1) = numel(V_AngList); nTriple(3,1) = numel(A_VgList); nTriple(4,1) = numel(A_VngList);
%%
% 
% BiasedCellIdx{1} = VgList;
% BiasedCellIdx{2} = VngList;
% BiasedCellIdx{3} = AgList;
% BiasedCellIdx{4} = AngList;

if OnlyInc == 1
    if OnlySelective == 1
        BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) == 1))];
        BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) == 1))];
        BiasedCellIdx{3} = [AgList(find(H_Data(3,AgList) == 1))];
        BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) == 1))];
    else
        BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) == 1)) VList(find(H_Data(1,VList) == 1))];
        BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) == 1)) VList(find(H_Data(2,VList) == 1))];
        BiasedCellIdx{3} = [AgList(find(H_Data(3,AgList) == 1)) AList(find(H_Data(3,AList) == 1))];
        BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) == 1)) AList(find(H_Data(4,AList) == 1))];
    end
elseif OnlyInc == 0
    BiasedCellIdx{1} = [VgList(find(H_Data(1,VgList) ~= 0)) VList(find(H_Data(1,VList) ~= 0))];
    BiasedCellIdx{2} = [VngList(find(H_Data(2,VngList) ~= 0)) VList(find(H_Data(2,VList) ~= 0))];
    BiasedCellIdx{3} = [AgList(find(H_Data(3,AgList) ~= 0)) AList(find(H_Data(3,AList) ~= 0))];
    BiasedCellIdx{4} = [AngList(find(H_Data(4,AngList) ~= 0))  AList(find(H_Data(4,AList) ~= 0))];
end

%%
C_Target = [5 6 5 6];
I_Target = [9 10; 7 8; 7 8; 9 10];

for i = 1:numel(BiasedCellIdx)
    clearvars Uni UniAmp
    Idx = BiasedCellIdx{i};
    for jj = 1:numel(Idx)
        Uni = nanmean(MergedSMIndAmp{i,Idx(jj)});
        Con = nanmean(MergedSMIndAmp{C_Target(i),Idx(jj)});
        Incon = nanmean([MergedSMIndAmp{I_Target(i,1),Idx(jj)}; MergedSMIndAmp{I_Target(i,2),Idx(jj)}]);
        
        % U_Data{2,i}(jj,1) = Uni;
        % C_Data{2,i}(jj,1) = Con;
        % I_Data{2,i}(jj,1) = Incon;
        if Abs == 1
            U_Data{2,i}(jj,1) = abs(Uni);
            C_Data{2,i}(jj,1) = abs(Con);
            I_Data{2,i}(jj,1) = abs(Incon);
        else
            U_Data{2,i}(jj,1) = (Uni);
            C_Data{2,i}(jj,1) = (Con);
            I_Data{2,i}(jj,1) = (Incon);
        end

        Uni = nanmean(MergedKMIndAmp{i,Idx(jj)});
        Con = nanmean(MergedKMIndAmp{C_Target(i),Idx(jj)});
        Incon = nanmean([MergedKMIndAmp{I_Target(i,1),Idx(jj)}; MergedKMIndAmp{I_Target(i,2),Idx(jj)}]);
        % 
        % U_Data{3,i}(jj,1) = Uni;
        % C_Data{3,i}(jj,1) = Con;
        % I_Data{3,i}(jj,1) = Incon;

        if Abs == 1
            U_Data{3,i}(jj,1) = abs(Uni);
            C_Data{3,i}(jj,1) = abs(Con);
            I_Data{3,i}(jj,1) = abs(Incon);
        else
            U_Data{3,i}(jj,1) = (Uni);
            C_Data{3,i}(jj,1) = (Con);
            I_Data{3,i}(jj,1) = (Incon);
        end
    end
end

%

for i = 1:3
    for j = 1:4
        U_Avg(i,j) = nanmean(U_Data{i,j}); U_Sem(i,j) = nanstd(U_Data{i,j})./sqrt(numel(U_Data{i,j}));
        C_Avg(i,j) = nanmean(C_Data{i,j}); C_Sem(i,j) = nanstd(C_Data{i,j})./sqrt(numel(C_Data{i,j}));
        I_Avg(i,j) = nanmean(I_Data{i,j}); I_Sem(i,j) = nanstd(I_Data{i,j})./sqrt(numel(I_Data{i,j}));
    end
end


%%

close all;

clearvars Color
Color = [197 197 197; 138 138 138; 0 0 0]./255;
MarkerColor = [255 0 0; 243 152 0; 0 0 255; 0 152 243]./255;


fig = figure('Position',[50 50 250 110]);
hold on

for j = 1:4
    for i = 1:3
        bar(j-0.25+(i-1)*0.25,C_Avg(i,j),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(j-0.25+(i-1)*0.25,C_Avg(i,j),C_Sem(i,j),C_Sem(i,j),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

for j = 1:4
    for i = 1:3
        scatter(j-0.25+(i-1)*0.25+0.05,U_Avg(i,j),10,'lineWidth',0.75,'MarkerEdgeColor',MarkerColor(j,:),'MarkerFaceColor','none');
        alpha(0.5);
        errorbar(j-0.25+(i-1)*0.25+0.05,U_Avg(i,j),U_Sem(i,j),U_Sem(i,j),'CapSize',0,'color',MarkerColor(j,:),'lineWidth',0.75)
        alpha(0.5);
    end
end

for j = 1:4
    for i = 1:3
        bar(4+j-0.25+(i-1)*0.25,I_Avg(i,j),0.25,'lineWidth',1,'EdgeColor',Color(i,:),'FaceColor','none');
        errorbar(4+j-0.25+(i-1)*0.25,I_Avg(i,j),I_Sem(i,j),I_Sem(i,j),'CapSize',3,'color',Color(i,:),'lineWidth',0.75)
    end
end

for j = 1:4
    for i = 1:3
        scatter(4+j-0.25+(i-1)*0.25+0.05,U_Avg(i,j),10,'lineWidth',0.75,'MarkerEdgeColor',MarkerColor(j,:),'MarkerFaceColor','none');
        alpha(0.5);
        errorbar(4+j-0.25+(i-1)*0.25+0.05,U_Avg(i,j),U_Sem(i,j),U_Sem(i,j),'CapSize',0,'color',MarkerColor(j,:),'lineWidth',0.75)
        alpha(0.5);
    end
end
xlim([0.4 8.6]); ylim([0 3]);
xticks([1:8]);
xticklabels({'VG','VNG','AG','ANG','VG','VNG','AG','ANG'});
set(gca,'TickDir','out','FontName','Arial','FontSize',6)
ylabel('Response amplitude, ΔF/F(%)','FontName','Arial','Fontsize',6);
xlabel('Type','FontName','Arial','Fontsize',6);

mkdir('Figure'); cd('Figure');
saveas(fig,'Fig 3h, Response Amp Change.svg');
cd ../

%%
