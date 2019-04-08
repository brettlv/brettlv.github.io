#!/usr/bin/python
# -*- coding: UTF-8 -*-

import re
import os
import xlrd,xlwt
import numpy as np
def main():
    x, y = [], []
    f = open('po_flux.dat','r')
    for line in f:
        line = line.rstrip()
        line = re.sub('\s+', ' ', line)
        line = line.lstrip(' ')
        line = line.rstrip(' ')
        line = line.split(' ')
        x.append(line[0])
        y.append(line[1])
    f.close()
    #return x,y
    #print(x[1], y[1])
    #readbook = xlrd.open_workbook('data.xls')
    writebook = xlwt.Workbook()
    sheet = writebook.add_sheet('1')
    #sheet = readbook.sheet_by_index(1)
    nrows = len(x)#行
    ncols = 2#列
    for i in np.arange(nrows):
        sheet.write(i,0,x[i])#写入excel，i行0列
        sheet.write(i,1,y[i])
    writebook.save('data.xls')



if __name__ == "__main__":
    main()
