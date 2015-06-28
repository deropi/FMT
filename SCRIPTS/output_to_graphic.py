# -*- coding: utf-8 -*-
"""
Created on Tue Jun 09 17:56:01 2015

@author: D
"""
import re
import csv
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
     AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.pyplot as plt

from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
     AnnotationBbox
from matplotlib.cbook import get_sample_data
from matplotlib.pyplot import figure, show

import numpy as np
import matplotlib.lines as mlines
from matplotlib.pyplot import xlim

omegas=[]
positions=[]
selection=[]
breakpoints=[]
file= open('output.txt','r')
content= file.readlines()
match= re.search('Significant positive selection sites found', content[1])
if match:
    positives=[]
    file= open('output.txt','r')
    content= file.readlines()
    content1= content[4:len(content)]
    reader=csv.reader(content1,delimiter='\t')
    for Position,Omega,Significant, bp in reader:
        print Position
        positions.append(Position)
        omegas.append(Omega)
        if bp== '*':
            breakpoints.append(Position)
        if float(Omega) >= 1:
           positives.append(Position) 
        if Significant == '*':
            selection.append(Position)
    fasta= open('protSEQ.fas', 'r')
    fasta1= fasta.readlines()
    aas= ''
    for i in range(1, len(fasta1)):
        match= re.search('>', fasta1[i])
        if match:
            break
        else:
            aas= aas + fasta1[i]
    frequencies= open('Positive_frequencies_h14_cluster.txt', 'w')
    frequencies.write('Pos\taa\tfreq')
    frequencies.write('\n')
    for i in range(0, len(positives)):
            posicion = int(positives[i])
            file = open('protSEQ.fas', 'r')
            aminoacidos=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q','R','S','T','U', 'V', 'W', 'Y']
            cuenta=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for line in file:
                if line[0] != '>':
                    for i in range(0, len(aminoacidos)):
                        if line[posicion -1] == aminoacidos[i]:
                            cuenta[i] = cuenta[i] +1
            count=0
            for x in range(0, len(cuenta)):
                count = count + cuenta[x]
                ultima=[]
            for y in range(0, len(cuenta)):
                if cuenta[y] != 0:
                    ultima.append(float(cuenta[y])/float(count))
                else:
                    ultima.append(0)
            for p in range(0, len(aminoacidos)):
                if ultima[p] != 0:
                    frequencies.write('%s\t%s\t(%s)' % (posicion, aminoacidos[p],ultima[p]))
                    frequencies.write('\n')
    frequencies.write('Averaged omega value\t')
    suma=0
    for i in range (0, len(omegas)):
        suma= suma + float(omegas[i])
    media= float(suma/len(omegas))
    frequencies.write('%s' %media)
    frequencies.write('\n')
    frequencies.write('Averaged omega value per partition\t')
    for i in range(0, len(breakpoints)):
        suma1=0
        if i ==0:
            for x in range(0, int(breakpoints[i])):
                suma1= suma1 + float(omegas[x])
            frequencies.write('%s' %(suma1/float(len(omegas[0:int(breakpoints[i])]))))
            frequencies.write('\t')
        else:
            for x in range(int(breakpoints[i-1]), int(breakpoints[i])):
                suma1= suma1 + float(omegas[x])
            frequencies.write('%s' %(suma1/float(len(omegas[int(breakpoints[i-1]):int(breakpoints[i])]))))
            frequencies.write('\t')
    
    frequencies.close()        
    if 1:
        fig, ax = plt.subplots()
        ax.plot(positions, omegas, linewidth=0.7)
        plt.xlim(int(positions[0]),int(len(positions)-1)+2)
        plt.ylim((float(min(omegas))-0.07), (float(max(omegas)))+0.2)
   
        ax.set_xlabel('Protein position')
        ax.set_ylabel('Omega value')
        for v in range(0, len(selection)):
                          ax.text((selection[v]), float(min(omegas)), '*', style='normal', color= 'red',  fontsize=7)
        for g in range(0, len(positions)):
            if float(omegas[g]) >1.05 and float(omegas[g]) < 1.3:
                ax.text((positions[g]), (float(min(omegas))-0.08)+0.03, aas[g], style='normal', color= 'lightsalmon',  fontsize=7)
            else:
                if float(omegas[g]) >= 1.3:
                    ax.text((positions[g]), (float(min(omegas))-0.08)+0.03, aas[g], style='normal', color= 'crimson',  fontsize=7)
                else:
                    if float(omegas[g]) >0 and float(omegas[g]) <= 0.5:
                        ax.text((positions[g]),(float(min(omegas))-0.08)+0.03, aas[g], style='normal', color= 'forestgreen',  fontsize=7)
                    else:
                        if float(omegas[g]) >=0.95 and float(omegas[g]) <= 1.03:
                            ax.text((positions[g]), (float(min(omegas))-0.08)+0.03, aas[g], style='normal',color= 'darkturquoise',  fontsize=7)
                        else:
                            ax.text((positions[g]), (float(min(omegas))-0.08)+0.03, aas[g], style='normal', color= 'lawngreen',  fontsize=7)
        matplotlib.pyplot.axhline(y=1,color= 'r', linewidth=0.4, linestyle="dashed")
        import matplotlib.patches as mpatches
        red_patch = mpatches.Patch(color='red', label='* Adaptive Positions')
        red_patch1= mpatches.Patch(color='crimson', label='Omega values >1.3')
        red_patch2= mpatches.Patch(color='lightsalmon', label='Omega values 1.05-1.3')
        red_patch3= mpatches.Patch(color='forestgreen', label='Omega values 0-0.5')
        red_patch4= mpatches.Patch(color='lawngreen', label='Omega values 0.5-0.95')
        red_patch5= mpatches.Patch(color='darkturquoise', label='Omega values 0.95-1.05 (neutral positions)')
        plt.legend(handles=[red_patch, red_patch1, red_patch2, red_patch3, red_patch4, red_patch5], loc= 2, prop={'size':7})
#        ax.annotate('* Positively selected positions', xy=(2, 1), xytext=(20, 1.72),color= 'red', fontsize= 10)
#        ax.annotate('Omega values >1.3', xy=(2, 1), xytext=(20, 1.68),color= 'crimson', fontsize= 10)
#        ax.annotate('Omega values 1.05-1.3', xy=(2, 1), xytext=(20, 1.64),color= 'lightsalmon', fontsize= 10)
#        ax.annotate('Omega values 0-0.5', xy=(2, 1), xytext=(20, 1.60),color= 'forestgreen', fontsize= 10)
#        ax.annotate('Omega values 0.5-0.95', xy=(2, 1), xytext=(20, 1.56),color= 'lawngreen', fontsize= 10)
#        ax.annotate('Omega values 0.95-1.05 (neutral positions)', xy=(2, 1), xytext=(20, 1.52),color= 'darkturquoise', fontsize= 10)


        plt.draw()
          
#plt.savefig('graphic.png', dpi=3000)
#plt.savefig('graphic.pdf', dpi=3000)
        plt.show()    
else:
    omegas=[]
    positions=[]
    file= open('output.txt','r')
    content= file.readlines()
    content1= content[4:len(content)]
    reader=csv.reader(content1,delimiter='\t')
    for Position,Omega, bp in reader:
        positions.append(Position)
        omegas.append(Omega)
        if bp == '*':
            breakpoints.append(Position)

    fasta= open('protSEQ.fas', 'r')
    fasta1= fasta.readlines()
    aas= ''
    for i in range(1, len(fasta1)):
        match= re.search('>', fasta1[i])
        if match:
            break
        else:
            aas= aas + fasta1[i]
            
    if 1:
        fig, ax = plt.subplots()
        ax.plot(positions, omegas, linewidth=0.7)
        plt.xlim(int(positions[0]),int(len(positions)-1)+2)
        plt.ylim((float(min(omegas))-0.015), (float(max(omegas)))+0.2)
        ax.set_xlabel('Protein position')
        ax.set_ylabel('Omega value')
        for g in range(0, len(positions)):
            if float(omegas[g]) >1.05 and float(omegas[g]) < 1.3:
                ax.text((positions[g]), (float(min(omegas))-0.015)+0.004, aas[g], style='normal', color= 'lightsalmon',  fontsize=7)
            else:
                if float(omegas[g]) >= 1.3:
                    ax.text((positions[g]), (float(min(omegas))-0.015)+0.004, aas[g], style='normal', color= 'crimson',  fontsize=7)
                else:
                    if float(omegas[g]) >0 and float(omegas[g]) <= 0.5:
                        ax.text((positions[g]), (float(min(omegas))-0.015)+0.004, aas[g], style='normal', color= 'forestgreen',  fontsize=7)
                    else:
                        if float(omegas[g]) >=0.95 and float(omegas[g]) <= 1.05:
                            ax.text((positions[g]), (float(min(omegas))-0.015)+0.004, aas[g], style='normal', color= 'darkturquoise',  fontsize=7)
                        else:
                            ax.text((positions[g]),(float(min(omegas))-0.015)+0.004, aas[g], style='normal', color= 'lawngreen',  fontsize=7)
        matplotlib.pyplot.axhline(y=1,color= 'r', linewidth=0.4, linestyle="dashed")
        import matplotlib.patches as mpatches
        red_patch1= mpatches.Patch(color='crimson', label='Omega values >1.3')
        red_patch2= mpatches.Patch(color='lightsalmon', label='Omega values 1.05-1.3')
        red_patch3= mpatches.Patch(color='forestgreen', label='Omega values 0-0.5')
        red_patch4= mpatches.Patch(color='lawngreen', label='Omega values 0.5-0.95')
        red_patch5= mpatches.Patch(color='darkturquoise', label='Omega values 0.95-1.05 (neutral positions)')
        plt.legend(handles=[red_patch1, red_patch2, red_patch3, red_patch4, red_patch5], loc= 2, prop={'size':7})
        frequencies= open('Positive_frequencies_h15_cluster.txt', 'w')
    frequencies.write('Averaged omega value\t')
    suma=0
    for i in range (0, len(omegas)):
        suma= suma + float(omegas[i])
    media= float(suma/len(omegas))
    frequencies.write('%s' %media)
    frequencies.write('\n')
    frequencies.write('Averaged omega value per partition\t')
    for i in range(0, len(breakpoints)):
        suma1=0
        if i ==0:
            for x in range(0, int(breakpoints[i])):
                suma1= suma1 + float(omegas[x])
            frequencies.write('%s' %(suma1/float(len(omegas[0:int(breakpoints[i])]))))
            frequencies.write('\t')
        else:
            for x in range(int(breakpoints[i-1]), int(breakpoints[i])):
                suma1= suma1 + float(omegas[x])
            frequencies.write('%s' %(suma1/float(len(omegas[int(breakpoints[i-1]):int(breakpoints[i])]))))
            frequencies.write('\t')
    
    frequencies.close()   


        #plt.draw()
        #fig.savefig('histoneh1_2_manual_result.pdf', format='pdf', dpi=1200)
plt.show()

fig.savefig('histoneh1_2_manual_result.pdf', format='pdf', dpi=1200)


