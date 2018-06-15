clear all
%%
BinIgnore=4;
BinFPC=dlmread('D:\171123_PZT calibration related\Some Result\20171123170755_128x_CurveEX_133-400_Bin_FPC.txt');

plot(BinFPC(:,1),BinFPC(:,2));
BinFPC_forAnalysis=BinFPC(BinIgnore:26,2);

PtPVariance=(max(BinFPC_forAnalysis)-min(BinFPC_forAnalysis))/mean(BinFPC_forAnalysis)
