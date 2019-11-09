function [dates,data,R_m,r_f] = loadStockData2(path)

data = xlsread(path,2,'B10:K1088');

date = xlsread(path,2,'A10:A1088');
date = num2str(date);
dates = datenum(date,'yyyymm');

r_f = xlsread(path,4,'C9:C1087');
R_m = xlsread(path,4,'B9:B1087');
