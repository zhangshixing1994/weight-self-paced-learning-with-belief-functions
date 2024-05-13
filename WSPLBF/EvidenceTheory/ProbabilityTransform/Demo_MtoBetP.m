clear all;
close all;
clc;

addpath(genpath(pwd));

[Data,label,Class_num] = NewDataSet(1);

mass= ECMBBA(Data,label,Data);
% mass=ECM(Data,3);

[PrPl] = MtoPrPl(mass,Class_num);
 [PrBl] = MtoPrBl(mass,Class_num);
[PraPl] = MtoPraPl(mass,Class_num);
 [BetP] = MtoBetP(mass);

prop=[0,0,0
    1,0,0
    0,1,0
    1/2,1/2,0
    0,0,1
    1/2,0,1/2
    0,1/2,1/2
    1/3,1/3,1/3];
 [MTP] = MyTransformP(mass,prop);