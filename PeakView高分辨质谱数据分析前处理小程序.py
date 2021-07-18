#Python3
#PeakView高分辨质谱数据分析前处理小程序

import numpy as np
import pandas as pd

Element = {'H':[[1.0078250322, 2.01410177811], [0.999885, 0.000115]],
               'Li':[[7.01600344, 6.015122887], [0.9241, 0.0759]],
               'B':[[11.00930517, 10.0129369], [0.801, 0.199]],
               'C':[[12.00000000,13.003354835], [0.9893, 0.0107]],
               'N':[[14.003074004, 15.000108899], [0.99636, 0.00364]],
               'O':[[15.994914619, 16.999131757, 17.999159613], [0.99757, 0.00038, 0.00205]],
               'F':[[18.998403163],[1]],
               'Na':[[22.98976928],[1]],
               'Mg':[[23.98504170, 25.9825930, 24.9858370], [0.7899, 0.1101, 0.1000]],
               'Al':[[26.9815384],[1]],
               'Si':[[27.976926535, 28.976494665, 29.9737701], [0.92223, 0.04685, 0.03092]],
               'P':[[30.973761999],[1]],
               'S':[[31.972071174, 33.9678670, 32.971458910, 35.967081], [0.9499, 0.0425, 0.0075, 0.0001]],
               'Cl':[[34.9688527, 36.9659026], [0.7576, 0.2424]],
               'K':[[38.96370649, 40.96182526, 39.9639982],[0.932581, 0.067302, 0.000117]],
               'Ca':[[39.9625909, 43.955482, 41.958618, 47.9525229, 42.958766, 45.95369], [0.96941, 0.02086, 0.00647, 0.00187, 0.00135, 0.00004]],
               'Br':[[78.918338, 80.916288], [0.5069, 0.4931]],
               'Ag':[[106.90509, 108.90476], [0.51839, 0.48161]],
               'I':[[126.90447],[1]]}
Electron = 0.000548579

def formula(fml):
    upchar = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    lowchar = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    num = ['0','1','2','3','4','5','6','7','8','9']
    position = []
    atomname = []
    atomnum = []
    for i in range(len(fml)):
        if fml[i] in upchar:
            position.append(i)
    for i in range(len(position)):
        if  position[i]+1 == len(fml) or fml[position[i]+1] in upchar:
            atomname.append(fml[position[i]])
            atomnum.append(1)
        elif fml[position[i]+1] in num:
            atomname.append(fml[position[i]])
            if position[i] == position[-1]:
                atomnum.append(int(fml[position[i]+1:]))
            #elif position[i]+1 == position[i+1]:
            #    atomnum.append(int(fml[position[i]+1]))
            #elif position[i]+1 < position[i+1]:
            else:
                atomnum.append(int(fml[position[i]+1:position[i+1]]))
        elif fml[position[i]+1] in lowchar:
            atomname.append(fml[position[i]:position[i]+2])
            if fml[position[i]+2] in num:
                atomnum.append(int(fml[position[i]+2:position[i+1]]))
            elif fml[position[i]+2] in upchar:
                atomnum.append(1)
            else:
                print("输入分子式错误 %s " %(fml))
                return None
    constitute = [atomname, atomnum]        
    return constitute

def exactms(constitute):
    ms = 0.0000
    for atomname, atomnum in zip(constitute[0], constitute[1]):
        ms = ms + Element[atomname][0][0]*atomnum
    return ms

def molwt(constitute):
    ms = 0.0000
    for atomname, atomnum in zip(constitute[0], constitute[1]):
        atomaverms = sum([w*m for w, m in zip(Element[atomname][0], Element[atomname][1])])
        ms = ms + atomaverms*atomnum
    return ms

def process():
    df1 = pd.read_excel(r"C:\Users\fangdoubleqiang\Desktop\HRMS.xlsx", sheet_name='InputInfo')    
    dfrow = len(df1)
    for i in range(dfrow):
        df1.loc[i,'No.'] = i+1
        df1.loc[i,'MoleculaWeight'] = round(molwt(formula(df1['MoleculaFormula'][i])), 2)
    print(df1)
    
    df2 = pd.read_excel(r"C:\Users\fangdoubleqiang\Desktop\HRMS.xlsx", sheet_name='OutputInfo')
    df2['No.'] = df1['No.']
    for i in range(dfrow):
        df2.loc[i, '[M+H]+'] = round(exactms(formula(df1['MoleculaFormula'][i])) + Element['H'][0][0] - Electron, 4)
        df2.loc[i, 'width+H'] = round(df2['[M+H]+'][i]*0.000015, 5)
        df2.loc[i, 'Formula+H'] = df1['MoleculaFormula'][i]
        df2.loc[i, '[M+Na]+'] = round(exactms(formula(df1['MoleculaFormula'][i])) + Element['Na'][0][0] - Electron, 4)
        df2.loc[i, 'width+Na'] = round(df2['[M+Na]+'][i]*0.000015, 5)
        df2.loc[i, 'Formula+Na'] = df1['MoleculaFormula'][i]
        df2.loc[i, '[M-H]-'] = round(exactms(formula(df1['MoleculaFormula'][i])) - Element['H'][0][0] + Electron, 4)
        df2.loc[i, 'width-H'] = round(df2['[M-H]-'][i]*0.000015, 5)
        df2.loc[i, 'Formula-H'] = df1['MoleculaFormula'][i]
        df2.loc[i, '[M+Cl]-'] = round(exactms(formula(df1['MoleculaFormula'][i])) + Element['Cl'][0][0] + Electron, 4)
        df2.loc[i, 'width+Cl'] = round(df2['[M+Cl]-'][i]*0.000015, 5)
        df2.loc[i, 'Formula+Cl'] = df1['MoleculaFormula'][i]
        
    print(df2)
    with pd.ExcelWriter(r"C:\Users\fangdoubleqiang\Desktop\HRMS.xlsx") as writer:
        df1.to_excel(writer, sheet_name='InputInfo', index=False)
        df2.to_excel(writer, sheet_name='OutputInfo', index=False)

process()

