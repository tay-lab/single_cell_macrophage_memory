%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure S3 related to
% Encoding innate immune memory in single macrophages through cell intrinsic and paracrine feedback
% Wang AG, Mercado Vasquez JG, Son M, et al. Corresponding author: tays@uchicago.edu
%
% Performance of the LSTM model and of a weighted random classifier (Figure S3A-C: tolerized LPS response,
% clusters of Figure 2; Figure S3D-F: primed TNF response, clusters of Figure 3).
% Uses the saved 5-fold cross-validation results in classifier_results/ (written by Fig2.m and Fig3.m; set
% runTraining = true there to retrain). A/D: ROC curves of the LSTM model (shaded: confidence interval across
% folds); B/E: F1 score of the random classifier per fold; C/F: its confusion matrix (columns normalized).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath functions simulation_results classifier_results data/flow   % run from the main code folder: all paths are relative to the current folder
files={"Fig2_classifiers.mat","Fig3_classifiers.mat"}; lab=["Figure S3A-C (CpG -> LPS)","Figure S3D-F (polyI:C -> TNF)"];
classNames=[{'1'} {'2'} {'3'} {'4'} {'5'}]; k=5;
for f=1:2
        R=load(files{f});   % found on the path (classifier_results/)
    cvLabels=cell(k,1); cvScores=cell(k,1);
    for fold=1:k
        cvLabels{fold}=R.datP{fold}{3};
        cvScores{fold}=R.datP{fold}{4};
    end
    rocObj=rocmetrics(cvLabels,cvScores,classNames);
    figure("Name",lab(f),"Position",[100 100 1200 350])
    subplot(1,3,1)                                   % A / D: LSTM ROC curves
    plot(rocObj,ShowConfidenceIntervals=true,ShowModelOperatingPoint=false)
    title("LSTM Model")
    subplot(1,3,2)                                   % B / E: random classifier F1 score
    boxplot(R.f1Random','Widths',0.5); hold on
    scatter(1:k,R.f1Random,300,".","blue")
    ylabel("F1 score"); xlabel("Group"); ylim([0 .7]); title("Random Classifier")
    subplot(1,3,3)                                   % C / F: random classifier confusion matrix
    totalConfMatrix=sum(R.CsRandom,3);
    cm=totalConfMatrix./sum(totalConfMatrix);
    h=heatmap(round(100*cm,1)); h.CellLabelFormat="%.1f%%"; h.ColorbarVisible="off";
    h.Title="Random Classifier"; h.XLabel="Predicted Group"; h.YLabel="True Group"; colormap(gca,sky)
end