function [dates,data,R_m,r_f] = loadStockData1(path)

data = xlsread(path,1,'B4:AE1088');

date = xlsread(path,1,'A4:A1088');
date = num2str(date);
dates = datenum(date,'yyyymm');

r_f = xlsread(path,4,'C3:C1087');
R_m = xlsread(path,4,'B3:B1087');
