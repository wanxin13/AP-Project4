function [dates,data,R_m,r_f] = loadStockData3(path)

data = xlsread(path,3,'B4:Z1088');

date = xlsread(path,3,'A4:A1088');
date = num2str(date);
dates = datenum(date,'yyyymm');

r_f = xlsread(path,4,'C3:C1087');
R_m = xlsread(path,4,'B3:B1087');
