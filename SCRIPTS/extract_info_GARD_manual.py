# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 20:28:02 2015

@author: D
"""

import re

positions=[]
aas=[]
omega1=[]
categoria=[]
positives=[]
breakpoints=[]
for exis in range(1, 3):
    
    file= open ("rst%i" %exis, "r")
    lnL= []
    archivo=[]
    positively_positions= []
    for line in file:
            if line.startswith("lnL"):
                lnL.append(line)   
            archivo.append(line)
    lnL2= [[],[]]
    for z in reversed(range(0, len(archivo))):
        match3 = re.search('Positively', archivo[z])
        if match3:
            for o in range (z, len(archivo)):
                positively_positions.append(archivo[o])
            break
    for i in range (0, len(lnL)):
            for c in range (0, len(lnL[i])):
                if lnL[i][c] == '0' or lnL[i][c] == '.' or lnL[i][c] == '1' or lnL[i][c] == '2' or lnL[i][c] == '3' or lnL[i][c] == '4' or lnL[i][c] == '5' or lnL[i][c] == '6' or lnL[i][c] == '7' or lnL[i][c] == '8' or lnL[i][c] == '9' or lnL[i][c] == '-':
                        lnL2[i].append(lnL[i][c])
    values=[]
    for x in range (0, len(lnL2)):
        t = ''.join(lnL2[x])
        values.append(float(t))
        
    LTR = (values[1] - values[0])*2
    omegas= []

    if LTR >= 4:
        for i in range (0, len(archivo)):
                if archivo[i] == lnL[1]:
                        for x in range (i, len(archivo)):
                                print archivo[x]
                                omegas.append(archivo[x])
                                ##print (archivo[x])
                                ##outfile.write(archivo[x])
        for z in range (0, len(omegas)):
                for a in range (0, len (omegas[z])):
                        match = re.search( '[A-Z|-]', omegas[z][a])
                        if match:
                                if omegas[z][10] != ' ' and omegas [z][0] == ' ' and a != 110:
                                        aas.append(omegas[z][a])
                                        #positions.append(omegas[z][1:4])
                                        omega1.append(omegas[z][103:108])
                                        categoria.append(omegas[z][99])
        positive=[]
        positive_site=[]
        breakpoints.append(len(omega1))
        for n in range (0, len(positively_positions)):
                        match4 = re.search( '[A-Z|-]', positively_positions[n])
                        if match4:
                            if positively_positions[n][8] == ' ':
                                print ('yes, it does')
                                positive.append(positively_positions[n][2:6])
                                positive_site.append(positively_positions[n][7])
        if exis-1 != 0:
            for i in range(0, len(positive)):
                positives.append(int(positive[i])+(breakpoints[exis-2]))
        else:
            for i in range(0, len(positive)): 
                positives.append(float(positive[i]))
    else:
        for i in range (0, len(archivo)):
                if archivo[i] == lnL[1]:
                        for x in range (i, len(archivo)):
                                omegas.append(archivo[x])
                                ##print (archivo[x])
                                ##outfile.write(archivo[x])
        for z in range (0, len(omegas)):
                for a in range (0, len (omegas[z])):
                        match = re.search( '[A-Z|-]', omegas[z][a])
                        if match:
                                if omegas[z][10] != ' ' and omegas [z][0] == ' ' and a != 110:
                                        aas.append(omegas[z][a])
                                        #positions.append(omegas[z][1:4])
                                        omega1.append(omegas[z][103:108])
                                        categoria.append(omegas[z][99])
        
        breakpoints.append(len(omega1))


for i in range(1, len(omega1)+1):
    positions.append('%i' %i)

if positives:
            out= open('results_h14_manual.txt', 'w')
            out.write('+-----------------------------------------------------+')
            out.write('\n')
            out.write('|      Significant positive selection sites found     |')
            out.write('\n')
            out.write('+-----------------------------------------------------+')
            out.write('\n')
            out.write('Position')
            out.write('\t')
            out.write('Omega')
            out.write('\t')
            out.write('Significant')
            out.write('\t')
            out.write('Breakpoint')
            out.write('\n')
            for i in range(0, len(positions)):
                alarma=0
                alarma2=0
                out.write('%s' %(positions[i]))
                out.write('\t')
                out.write('%s' %(omega1[i])) 
                out.write('\t')
                for x in range(0, len(positives)):
                    if positives[x] == int(positions[i]):
                        out.write('*')
                        out.write('\t')
                        alarma=1
                if alarma==0:
                    out.write('-')
                    out.write('\t')
                for x in range(0, len(breakpoints)):
                    if breakpoints[x] == int(positions[i]):
                        out.write('*')
                        alarma2=1
                if alarma2==0:
                    out.write('-')                  
                out.write('\n')
            
            out.close()
elif not positives:
            out= open('results_h11_manual.txt', 'w')            
            out.write('+-----------------------------------------------------+')
            out.write('\n')
            out.write('|             No Positive selection found             |')
            out.write('\n')
            out.write('+-----------------------------------------------------+')
            out.write('\n')
            out.write('Position')
            out.write('\t')
            out.write('Omega')
            out.write('\t')
            out.write('Breakpoint')
            out.write('\n')
            for i in range(0, len(positions)):
                out.write('%s' %(positions[i]))
                out.write('\t')
                out.write('%s' %(omega1[i])) 
                out.write('\t')
                for x in range(0, len(breakpoints)):
                    if breakpoints[x] == positions[i]:
                        out.write('*')
                        break
                else:
                        alarma=0
                if alarma==0:
                    out.write('-')
                out.write('\n')

            out.close()
