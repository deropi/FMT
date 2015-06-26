#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from sys import argv
import getopt
import urllib
import os
import subprocess
import re






def align():
    os.system('/home/devani/PIPE/programs/clustalw2 -INFILE=/home/devani/PIPE/%s/protSEQ.fas -TYPE=PROTEIN -OUTFILE=/home/devani/PIPE/%s/clustalw.fasta -OUTPUT=fasta -OUTORDER=INPUT > /home/devani/PIPE/test/mierda' %(fold[0], fold[0]))
    
def getSEQ():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:o", ["help", "input=", "output="])
    except getopt.GetoptError, err:
        print str(err)
        #usage()
        sys.exit(2)
    archivo=[]
    #directoryPath = 'test1'
    #os.system("mkdir test1")
    for o in args:        
        out= o
        os.system("mkdir %s" %out)
    
    for o, a in opts:
        if o in ("h","--help"):
            print ("""
      --help      Print this
      --input     Input file with Accession Numbers required
      --output    Output fold name required.
     """)
            break
        elif o in ("-i", "--input"):
            t = open(a, 'r')
            for line in t:
                f = urllib.urlopen('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=%d&amp;rettype=fasta' %int(line))
                result = f.read().decode('utf-8')
                archivo.append(result)
                #print(result)
                os.system("cd %s" %out)
                outs = open('/home/devani/PIPE/%s/secuencias.fas' %out, "w+")
                for i in range(0, len(archivo)):
                    outs.write(archivo[i])
            #import translate_dna
            #align()
    fold.append(out)
    
#Ahora mismo hemos colocado el fasta con las secuencias como output, pero despu?s el out ser? un report con los datos, no ?sto...habr? que introducirlo dentro de otra rutina...      
        #elif o in ("o", "--output"):
        #   print ('luego habr? un out')
fold=[]
getSEQ()


#AQUÍ HAY QUE METER ALGO QUE ALINEE DNA, Y EN BASE A ESTE ALINEAMIENTO ESCOJA EL ATG. TENGO secuencias.fas.
def align_dna_firststep():
    os.system('/home/devani/PIPE/programs/clustalw2 -INFILE=/home/devani/PIPE/%s/secuencias.fas -TYPE=DNA -OUTFILE=/home/devani/PIPE/%s/clustalw_dna.fasta -OUTPUT=fasta -OUTORDER=INPUT > /home/devani/PIPE/%s/mierda' %(fold[0], fold[0], fold[0]))
    dna= open('/home/devani/PIPE/%s/clustalw_dna.fasta' %fold[0], 'r')
    secuencia_alineada= open ('/home/devani/PIPE/%s/secuencia_alineada.fas' %fold[0], 'w')
    header1 = []
    sequence1=[]
    header=''
    sequence=''
    for line in dna:
        if line[0] == ">":
            if header != '':
                header1.append(line)
                sequence1.append(sequence)
            header = line.strip()
            sequence = ''
        else:
            sequence += line.strip()
    sequence1.append(sequence)
    header1.append(line)
    start=0
    comienzo=start
    inicio=0
#    start= sequence1[0].find('ATG')
    for i in range(0, len(sequence1[0]),3):
        start= sequence1[0][inicio:].find('ATG')
        alarma=0
        for x in range(1, len(sequence1)):
            if sequence1[x][int(start)+inicio:int(start)+inicio+3] =='ATG':
                alarma=0
            else:
                alarma=1
                inicio= inicio + int(start)+3
                break
                #i=start+3
        if alarma==0:
            comienzo=int(start)
            break
    inicio = inicio + int(start)
    for i in range(0, len(sequence1)):
        secuencia_alineada.write(header1[i])
        secuencia_alineada.write('\n')
        for x in range(inicio, len(sequence1[i])):
            if sequence1[i][x] != '-':
                secuencia_alineada.write(sequence1[i][x])
        secuencia_alineada.write('\n')

    secuencia_alineada.close()
    dna.close()
    
    
    
    
align_dna_firststep()   









protSEQ= open('/home/devani/PIPE/%s/protSEQ.fas' %fold[0], 'w')
dnaSEQ= open('/home/devani/PIPE/%s/dnaSEQ.fas' %fold[0], 'w')
ultima_prot=[]
ultima_dna=[]
def translate_dna(sequence):

    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    proteinsequence = ''
    dnasequence = ''
    #start = sequence.find('ATG')
    #sequencestart = sequence[int(start):]
    sequencestart = sequence[0:]
    cds = str(sequencestart[:len(sequence)])

    for n in range(0,len(cds),3):
        if cds[n:n+3] in codontable:
            if codontable[cds[n:n+3]] == '_' :
                break
            else:
                
                proteinsequence += codontable[cds[n:n+3]]
                dnasequence += cds[n:n+3]
        sequence = ''
    ultima_prot.append(proteinsequence)
    ultima_dna.append(dnasequence)
    protSEQ.write(proteinsequence)
    dnaSEQ.write(dnasequence)
    protSEQ.write('\n')
    dnaSEQ.write('\n')


inputsequ = open('/home/devani/PIPE/%s/secuencia_alineada.fas'%fold[0], 'r')
header = ''
sequence = ''
for line in inputsequ:
    if line[0] == ">":
        if header != '':
            protSEQ.write(header)
            protSEQ.write('\n')
            dnaSEQ.write(header)
            dnaSEQ.write('\n')
            translate_dna(sequence)

        header = line.strip()
        sequence = ''
    else:
        sequence += line.strip()
        
protSEQ.write(header)
protSEQ.write('\n')
dnaSEQ.write(header)
dnaSEQ.write('\n')
protSEQ.write(ultima_prot[len(ultima_prot)-1])
dnaSEQ.write(ultima_dna[len(ultima_dna)-1])

protSEQ.close()
dnaSEQ.close()
inputsequ.close()

align()

proteina_lista=[]
dna_lista=[]
#no funciona porque clustalw.fasta y dnaSEQ.fas no están en el mismo orden. Ya funciona!!!
proteina_alineada = open('/home/devani/PIPE/%s/clustalw.fasta' %fold[0], 'r')
dna_alineada = open('/home/devani/PIPE/%s/dnaSEQ.fas' %fold[0], 'r+')
    
header1 = ''
sequence1 = ''
for line in proteina_alineada:
    if line[0] == ">":
        if header1 != '':
            proteina_lista.append(sequence1)
        header1 = line.strip()
        sequence1 = ''
    else:
        sequence1 += line.strip()
            
proteina_lista.append(sequence1)
header2=''
sequence2=''
for line in dna_alineada:
    if line[0] == ">":
        if header2 != '':
            dna_lista.append(sequence2)
        header2 = line.strip()
        sequence2 = ''
    else:
        sequence2 += line.strip()
dna_lista.append(sequence2)

proteina_alineada.close()
dna_alineada.close()



def aligned_prot_to_dna():
    nueva_dna_lista=[]
    for i in range(0, len(proteina_lista)):
        nueva_dna=''
        count = 0
        nueva_dna_lista_2=[]
        for x in range(0, len(proteina_lista[i])):
            a = (x*3) - (count*3)
            if proteina_lista[i][x] == '-':
                count = count +1
                a = (x*3) - (count*3)
                nueva_dna = ('---')
            else:
                nueva_dna = dna_lista[i][a:(a+3)]
            nueva_dna_lista_2.append(nueva_dna)
        nueva_dna_lista.append(''.join(nueva_dna_lista_2))

    dna_aligned = open ('/home/devani/PIPE/%s/dna_clustalw.fasta' %fold[0], 'w')
    proteina_alineada = open('/home/devani/PIPE/%s/clustalw.fasta' %fold[0], 'r')
    count2= 0
    for line in proteina_alineada:

        if line[0] == '>':
            dna_aligned.write(line)
            dna_aligned.write(nueva_dna_lista[count2])
            dna_aligned.write('\n')
            count2 = count2 +1 
    dna_aligned.close()
    proteina_alineada.close()
    
        
    
#tienes que cerrar los files!!!!!!!

aligned_prot_to_dna()

def gaps():
    sin_gaps= open('/home/devani/PIPE/%s/dna_aligned_no_gaps.fas' %fold[0], 'w')
    con_gaps= open ('/home/devani/PIPE/%s/dna_clustalw.fasta'%fold[0], 'r')
    con_gaps_2= open ('/home/devani/PIPE/%s/dna_clustalw.fasta' %fold[0], 'r')
    gaps= [] 
    for line in con_gaps:
        if line[0] != '>':
            for i in range (0, len(line)):
                if line[i] == '-':
                    gaps.append(i)
    for line in con_gaps_2:
        if line[0] == '>':
            sin_gaps.write(line)
        else:           
            for k in range(0, len(line)):
                contador= 0
                for x in range(0, len(gaps)):
                    if k == gaps[x]:
                        contador = contador +1
                        break
                if contador == 0:
                    sin_gaps.write(line[k])

    sin_gaps.close()
    con_gaps.close()
    con_gaps_2.close()
                        

gaps()

#model_code= ''
def NJ_modeltest():
    settings_tree=open ('/home/devani/PIPE/%s/settings_arbol' %fold[0], 'w')
    settings_tree.write('1\n2\n1\n/home/devani/PIPE/%s/dna_aligned_no_gaps.fas\n1\n1\n1\ny\n/home/devani/PIPE/%s/arbol_inicial.nex\n' %(fold[0], fold[0]))
    settings_tree.close()
    settings_modeltest=open('/home/devani/PIPE/%s/modeltest.txt' %fold[0], 'w')
    settings_modeltest.write('/home/devani/PIPE/%s/dna_aligned_no_gaps.fas\n/home/devani/PIPE/%s/arbol_inicial.nex\n4\n2\n/home/devani/PIPE/%s/model_test_out\n' %(fold[0], fold[0], fold[0]))
    settings_modeltest.close()
    os.chdir('/home/devani/PIPE/programs')
#dándole una salida a un directorio que luego borraremos, evitamos que corra directamente en línea de comandos. conseguimos solo el out. revisar donde están el settings_arbol y el modeltest.txt y los parámetros
    os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/NeighborJoining.bf < /home/devani/PIPE/%s/settings_arbol > /home/devani/PIPE/%s/mierda' %(fold[0], fold[0]))
    cmd = ("sed -i 's/|//g' /home/devani/PIPE/%s/dna_aligned_no_gaps.fas" %fold[0])
    os.system(cmd)
    cmd= ("sed -i 's/[.]//g' /home/devani/PIPE/%s/dna_aligned_no_gaps.fas" %fold[0])
    os.system(cmd)
    cmd = ("sed -i 's/|//g' /home/devani/PIPE/%s/arbol_inicial.nex " %fold[0])
    os.system(cmd)
    cmd= ("sed -i 's/[.]//g' /home/devani/PIPE/%s/arbol_inicial.nex" %fold[0])
    os.system(cmd)
    os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/ModelTest.bf < /home/devani/PIPE/%s/modeltest.txt > /home/devani/PIPE/%s/out_run_jmodeltest.txt ' %(fold[0], fold[0]))
#meter este modeltest_out en el directorio basura???    
    models_out= open('/home/devani/PIPE/%s/out_run_jmodeltest.txt' %fold[0], 'r')
    model_code= ''
    for line in models_out:
        if line[0:13] == 'Model String:':
            model_code = line[13:20]
            break
    settings_GARD= open('/home/devani/PIPE/%s/settings_GARD.txt' %fold[0], 'w')
    settings_GARD.write('/home/devani/PIPE/%s/dna_aligned_no_gaps.fas\n%s1\n/home/devani/PIPE/%s/GARD.txt\n' %(fold[0], model_code, fold[0]))
    settings_GARD.close()


    
            
    
#model_code= ''   
NJ_modeltest()
#print model_code

def GARD():
#el file en el GTAcluster quetiene los parámetros del GARD es un settingsque se encuentra en el directorio que se genera cada vez
#cada vez que corremos el GARD salen diferentes resultados...es así...no la estás cagando tanto...
    os.system('mpiexec.hydra -n 4 /home/devani/PIPE/programs/HYPHYMPI /home/devani/PIPE/programs/TemplateBatchFiles/GARD.bf < /home/devani/PIPE/%s/settings_GARD.txt > /home/devani/PIPE/%s/results_GARD.txt' %(fold[0], fold[0]))
    #os.system('/home/devani/PIPE/programs/HYPHYMPI')
    #os.system('/home/devani/PIPE/programs/TemplateBatchFiles/GARD.bf < /home/devani/PIPE/test/dna_aligned_no_gaps.fas')
    GARD_out= open('/home/devani/PIPE/%s/GARD.txt_finalout' %fold[0], 'r')
    lines=GARD_out.readlines()
    #breakpoints=[]
    breakpoints_lines=[]
    for i in range(0, len(lines)):
        match= re.search('BEGIN ASSUMPTIONS;', lines[i])
        if match:
            breakpoints_lines.append(lines[i+1:len(lines)])
    for i in range(0, len(breakpoints_lines[0])):
        for x in range(0, len(breakpoints_lines[0][i])):
            
            if breakpoints_lines[0][i][x] == '-':
                breakpoints.append(breakpoints_lines[0][i][x+1:(len(breakpoints_lines[0][i])-2)])
breakpoints=[]    
GARD()



breakpoints1=[]
def Partitions():

    secuencia= open('/home/devani/PIPE/%s/dna_aligned_no_gaps.fas' %fold[0], 'r')
    dna_alineado=[]
    cabeceras=[]
    header=''
    sequence=''
    for line in secuencia:
        if line[0] == ">":
            cabeceras.append(line)
            if header != '':
                dna_alineado.append(sequence)
            header = line.strip()
            sequence = ''
        else:
            sequence += line.strip()
    dna_alineado.append(sequence)
    count=1
#Los puntos del GARD no coinciden con final de codones. breakpoints1 es una lista de los puntos en que corta
#las secuencias teniendo esto en cuenta

    #breakpoints1=[]
    if int(breakpoints[0])%3 == 0:
        breakpoints1.append(breakpoints[0])
    elif (int(breakpoints[0])-1)%3 == 0:
        breakpoints1.append(int(breakpoints[0])-1)
    elif (int(breakpoints[0])-2)%3 == 0:
        breakpoints1.append(int(breakpoints[0])-2)
    for i in range(1, len(breakpoints)):
        if ((int(breakpoints[i])-int(breakpoints1[i-1]))%3) == 0:
            breakpoints1.append(int(breakpoints[i]))
        elif (((int(breakpoints[i])+1)-int(breakpoints1[i-1]))%3) == 0:
            breakpoints1.append(int(breakpoints[i])+1)
        elif (((int(breakpoints[i])+2)-int(breakpoints1[i-1]))%3) == 0:
            breakpoints1.append(int(breakpoints[i])+2)           
    for i in range(0, len(breakpoints1)):
        particion= open('/home/devani/PIPE/%s/dna_aligned%i.txt' %(fold[0], count), 'w')
        if i == 0:
            for x in range(0, len(cabeceras)):
                particion.write(cabeceras[x])
                particion.write('\n')
                particion.write(dna_alineado[x][0:int(breakpoints1[0])+1])
                particion.write('\n')
#para que no corte el final de la secuencia
        elif i == (len(breakpoints)-1):
            for x in range(0, len(cabeceras)):
                particion.write(cabeceras[x])
                particion.write('\n')
                a= int(breakpoints1[i-1])
                b= len(dna_alineado[x])+1
                particion.write(dna_alineado[x][a:b])
                particion.write('\n')
                
        else:
            for x in range(0, len(cabeceras)):
                particion.write(cabeceras[x])
                particion.write('\n')
                a= int(breakpoints1[i-1])
                b= int(breakpoints1[i])+1
                particion.write(dna_alineado[x][a:b])
                particion.write('\n')

                        
                        
        particion.close()
        count= count+1
        fragmentos.append(count)

        

#frahmentos es una lista que contiene el numero de los fragmentos...sirve para buscar luego los archivos que se generaron
fragmentos = []
#hacer partitions solo si hay bp. Necesitaria ver como quedan los outs del GARD para esta condicion

if len(breakpoints) != 1:
    Partitions()
else:
    fragmentos=[]




#Positive Selection NielsenYang.bfUse a series of REL models (with constant ) to
#test for selection. This analysis (with the GY94
#model) is identical to PAML analyses, e.g. M8a
#v M8b tests.
 
#settings_arbol= open('/home/devani/PIPE/test/settings_arbol', 'w')


def NJ_Model_foreach_fragment():
    count=1
    count2=1
    for x in range(0, (len(fragmentos))):
        settings_arbol= open('/home/devani/PIPE/%s/settings_arbol' %fold[0], 'w')
        settings_arbol.write('1\n2\n1\n/home/devani/PIPE/%s/dna_aligned%i.txt\n1\n1\n1\ny\n/home/devani/PIPE/%s/arbol_primero%i\n' %(fold[0], count, fold[0],count2) )
        settings_arbol.close()
        os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/NeighborJoining.bf < /home/devani/PIPE/%s/settings_arbol > /home/devani/PIPE/%s/mierda' %(fold[0], fold[0]))
        selection_settings=open('/home/devani/PIPE/%s/settings_selection' %fold[0], 'w')
        #modelo 16 vs. 12: beta&1, beta&normal>1, hay que cambiar el cutoff a0.5. La posterior cutoff es el BEB de PAML. n a 12 porque en modelo &1 hay 11 aquí, no 10.Pero en principio, los parámetros no varían con el número de categorias.
        #modelos 8 vs.9 son beta y beta&w :es decir, M7 Y M8
        selection_settings.write('1\n/home/devani/PIPE/%s/dna_aligned%i.txt\n/home/devani/PIPE/%s/arbol_primero%i\n2\n8\n9\nd\n1\n3\n0.5\n10\n/home/devani/PIPE/%s/results_selection%i\n'%(fold[0], count, fold[0], count, fold[0], count))
        #selection_settings.write('1\n/home/devani/PIPE/test/dna_aligned%i.txt\n/home/devani/PIPE/test/arbol_primero%i\n2\n16\n13\nd\n1\n3\n0.5\n11\n/home/devani/PIPE/test/results_selection_segundaprueba%i\n'%(count, count, count))
        selection_settings.close()
        os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/NielsenYang.bf < /home/devani/PIPE/%s/settings_selection > /home/devani/PIPE/%s/mierda' %(fold[0], fold[0]))
        count = count +1
        count2= count2 +1

 

def NJ_Model_foreach_fragment_for_no_recombination():
        count=1
        settings_arbol= open('/home/devani/PIPE/%s/settings_arbol' %fold[0] , 'w')
        settings_arbol.write('1\n2\n1\n/home/devani/PIPE/%s/dna_aligned.txt\n1\n1\n1\ny\n/home/devani/PIPE/%s/arbol_primero\n'  %(fold[0], fold[0]))
        settings_arbol.close()
        os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/NeighborJoining.bf < /home/devani/PIPE/%s/settings_arbol > /home/devani/PIPE/%s/mierda' %(fold[0], fold[0]))
        selection_settings=open('/home/devani/PIPE/%s/settings_selection' %fold[0], 'w')
        #modelo 16 vs. 12: beta&1, beta&normal>1, hay que cambiar el cutoff a0.5. La posterior cutoff es el BEB de PAML. n a 12 porque en modelo &1 hay 11 aquí, no 10.Pero en principio, los parámetros no varían con el número de categorias.
        selection_settings.write('1\n/home/devani/PIPE/%s/dna_aligned.txt\n/home/devani/PIPE/%s/arbol_primero%i\n2\n16\n13\nd\n1\n3\n0.5\n10\n/home/devani/PIPE/%s/results_selection%i\n' %(fold[0], fold[0], count, fold[0], count))
        #selection_settings.write('1\n/home/devani/PIPE/test/dna_aligned%i.txt\n/home/devani/PIPE/test/arbol_primero%i\n2\n16\n13\nd\n1\n3\n0.5\n11\n/home/devani/PIPE/test/results_selection_segundaprueba%i\n'%(count, count, count))
        selection_settings.close()
        os.system('mpiexec.hydra -n 4 ./HYPHYMPI ./TemplateBatchFiles/NielsenYang.bf < /home/devani/PIPE/%s/settings_selection%i > /home/devani/PIPE/%s/mierda' %(fold[0], count, fold[0]))
    
if len(breakpoints) != 1:    
    NJ_Model_foreach_fragment()
else:
    NJ_Model_foreach_fragment_for_no_recombination()

#fuera de las subrutinas deberíamos poner el valor del contador que nos ha generado los archivos del GARD
omegas= []
positions=[]
#limites tiene el numero de posiciones de cada file 
limites=[]
limites_final=[]
#list_positives es una lista que contiene dentro listas con los sitos positivos. se han ido añadiendo en orden
#de output del HYPHY. ya se han sumado las posiciones anteriores. con lo cual está en orden
list_positives=[]
list_positives_final=[]
lnL_out=[]
def extract_info():
    count=1
    for i in range(1, len(breakpoints)+1):
        ln= []
        HYPHY = open('/home/devani/PIPE/%s/results_selection%s' %(fold[0], i), 'r')
        for line in HYPHY:
            if line.startswith('-') and line[1] != '-':
                ln.append(line)
        lnL= (float(ln[0][1:10])- float(ln[1][1:10]))*2
        lnL_out.append(float(lnL))
        HYPHY.close()
#CAMBIAR EL VALOR
#==============================================================================
# Aqupara coger los datos de posible seleccion positiva        
#==============================================================================
        if lnL > 5:
            rates=[]
            probabilities=[]
            count1 = 0
            f= open('/home/devani/PIPE/%s/results_selection%s' %(fold[0], i), 'r')
            content = f.readlines()
            f.close()
            for x in range (0, len(content)):
                match = re.search('(MODEL 8)', content[x])
                count1 = count1 +1
                if match:
                    break
            p= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'w')
            for x in range (count1-1, len(content)):
                match2 = re.search('MODEL 8', content[x])
                match10 = re.search('MODEL 7', content[u])
                if match2:
                    for u in range(x, len(content)):
                        match10 = re.search('MODEL 7', content[u])
                        if match10:
                            break
                        else:
                            p.write(content[u])
                elif match10:
                    break
            p.close()
            g= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'r')
            otro= open('/home/devani/PIPE/%s/otro.txt'  %fold[0], 'w')
            contenido = g.readlines()
            for x in range (0, len(contenido)):
                if len(contenido[x]) > 70:
                    otro.write(contenido[x])
            g.close()
            otro.close()
            cuenta3=0
            j= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'r')
            content = j.readlines()
            for x in range(0, len(content)):
                match3= re.search('Sites with dN/dS>1', content[x])
                if match3:
                    break
                else:
                    cuenta3= cuenta3+1
            j.close()
            positivos=[]            
            cuenta3= cuenta3+1
            modelos= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'r')
            contenido2= modelos.readlines()
            for x in range(cuenta3, len(contenido2)):
                match4= re.search('------------------------------------------------', contenido2[x])
                if match4:
                    break
                else:
                    positivos.append(contenido2[x])
            modelos.close()
            positives_extracted=[]
            for x in range(0, len(positivos)):
                extract= ''
                for u in range(0, len(positivos[x])):
                    if positivos[x][u] != ' '  and positivos[x][u] != '\n':
                        extract= extract + positivos[x][u]
                    else:
                        break
                if extract != '':
                    positives_extracted.append(''.join(extract))
            #for x in range(0, len(positives_extracted)):
            #    print positives_extracted[x]
#comprobar que pasa si el modelo es positivo pero no hay sitios +. creo que no coger?a nada...asi que funciona igual...no s?...comprobar
            list_positives.append(positives_extracted)                                
            file = open('/home/devani/PIPE/%s/otro.txt' %fold[0], 'r')
            for line in file:

                if line.startswith('R'):
                    rates.append(line.split('   '))
                else:
                    probabilities.append(line.split(' '))
            file.close()
            rates1= []
            rates2=''
            for x in range(0, len (rates[0][0])):
                match25= re.search('.', rates2)
                match= re.search('[0-9|.]', rates[0][0][x])
                if match:
                    rates2= rates2 + rates[0][0][x]
                elif rates[0][0][x] == '\t' or rates[0][0][x] == '\n':
                    if rates2 and match25:
                       rates1.append(rates2)
                    rates2=''
                elif rates[0][0][x] == '-':
                    rates2='0'
                    rates1.append(rates2)
                    rates2=''
                    

            match25= re.search('.', rates2)
            if rates2 != '' and match25:
                rates1.append(rates2)
            elif rates2 != '':
                rates2= '0'
                rates1.append(rates2)
            rates_final=[]
            for x in range(0, len(rates1)):
                if rates1[x] != '':
                    rates_final.append(float(rates1[x]))

            probabilities_final=[]
            for x in range(0, len(probabilities)):
                probs=''
                for u in range(0, len(probabilities[x][0])):
                    match1= re.search('[0-9|.]', probabilities[x][0][u])
                    match7= re.search('[e|-]', probabilities[x][0][u])
                    if match1 or match7:
                        probs = probs +  probabilities[x][0][u]
                    elif probabilities[x][0][u] == '\t' or probabilities[x][0][u] == '\n':
                        match8= re.search('e', probs)
                        if match8:
                            probs= '0'
                            probabilities_final.append(float(probs))
                            probs=''
                        else:    
                            probabilities_final.append(float(probs))
                            probs=''
                #probabilities_final.append(probs)
            probs1=[]
            for x in range(0, len(probabilities_final), len(rates_final)+1):
                probs=[]
                probs.append(probabilities_final[x:x+(len(rates_final)+1)])
                probs1.append(probs)
            averaged_omega=[]
            cuenta_posiciones_porfile=0
            for x in range(0, len(probs1)):                                                    
                media= 0 
                for u in range(1, len(probs1[0][0])):
                    media= media + float(probs1[x][0][u] * float(rates_final[u-1]))
                averaged_omega.append(media)
                cuenta_posiciones_porfile= cuenta_posiciones_porfile +1
            omegas.append(averaged_omega)
            limites.append(cuenta_posiciones_porfile)
            

            
            
#==============================================================================
#             aqu para no seleccion positiva
#==============================================================================
        elif lnL < 5:
            rates=[]
            probabilities=[]
            count1 = 0
            f= open('/home/devani/PIPE/%s/results_selection%s' %(fold[0], i), 'r')
            content = f.readlines()
            f.close()
            for x in range (0, len(content)):
                match = re.search('MODEL 7', content[x])
                count1 = count1 +1
                if match:
                    break
            p= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'w')
            for x in range (count1, len(content)):
                match2 = re.search('MODEL 8', content[x])
                if match2:
                    break
                else:
                    p.write(content[x])
            p.close()
            g= open('/home/devani/PIPE/%s/modelos.txt' %fold[0], 'r')
            otro= open('/home/devani/PIPE/%s/otro.txt' %fold[0], 'w')
            contenido = g.readlines()
            for x in range (0, len(contenido)):
                if len(contenido[x]) > 70:
                    otro.write(contenido[x])
            g.close()
            otro.close()

            
        
#hasta aqui tengo una lista... en la que tengo las lineas del file con los sitios positivos
#hay que sacar los valores y meterlos en una lista donde esten ordenados para la secuencia entera
            
                    
            file = open('/home/devani/PIPE/%s/otro.txt' %fold[0], 'r')
            for line in file:
                if line.startswith('R'):
                    rates.append(line.split('  '))
                else:
                    probabilities.append(line.split(' '))
            rates1= []
            rates2=''
            file.close()
            for x in range(0, len (rates[0][0])):
                match= re.search('[0-9|.]', rates[0][0][x])
                if match:
                    rates2= rates2 + rates[0][0][x]
                elif rates[0][0][x] == '\t' or rates[0][0][x] == '\n':
                    print rates2
                    if rates2 and float(rates2)<1:
                       rates1.append(rates2)
                    rates2=''
                elif rates[0][0][x] == '-':
                    rates2='0'
                    rates1.append(rates2)
                    rates2=''
            match25= re.search('.', rates2)
            if rates2 != '' and match25:
                rates1.append(rates2)
            elif rates2 != '':
                rates2= '0'
                rates1.append(rates2)
            rates_final=[]
            for x in range(0, len(rates1)):
                if rates1[x] != '':
                    rates_final.append(float(rates1[x]))

#                if rates2 != '':
#                    rates3.append(''.join(rates2)) 

            probabilities_final=[]
            for x in range(0, len(probabilities)):
                probs=''
                for u in range(0, len(probabilities[x][0])):
                    match7= re.search('[e|-]', probabilities[x][0][u])
                    match1= re.search('[0-9|.]', probabilities[x][0][u])
                    if match1 or match7:
                        probs = probs + (probabilities[x][0][u])
                    elif probabilities[x][0][u] == '\t' or probabilities[x][0][u] == '\n':
                        match8= re.search('e', probs)
                        if match8:
                            probs= '0'
                            probabilities_final.append(float(probs))
                            probs=''
                        else:    
                            probabilities_final.append(float(probs))
                            probs=''
                #probabilities_final.append(probs)
            probs1=[]
            averaged_omega=[]
            cuenta_posiciones_porfile=0
            for x in range(0, len(probabilities_final), len(rates_final)+1):
                probs=[]
                probs.append(probabilities_final[x:x+(len(rates_final)+1)])
                probs1.append(probs)
                averaged_omega=[]
                cuenta_posiciones_porfile=0
            for x in range(0, len(probs1)):
                media= 0 
                for u in range(1, len(probs1[0][0])):
                    media= media + (float(probs1[x][0][u] * float(rates_final[u-1])))
                averaged_omega.append(media)
                cuenta_posiciones_porfile= cuenta_posiciones_porfile +1
            omegas.append(averaged_omega)
            limites.append(cuenta_posiciones_porfile)
            list_positives.append('')
        count=count+1   

#    for i in range(0, len(omegas)):
#        for x in range(0, len(omegas[i])):
#            positions.append(count)
#            count = count+1
#    for i in range(1, len(list_positives)):
#        for x in range(0, len(list_positives[i])):
#            list_positives[i][x] = int(list_positives[i][x]) + limites[i-1]

        
        
            
extract_info()

count5=0
for i in range(0, len(omegas)):
    for x in range(0, len(omegas[i])):
        count5= count5+1

for i in range(1, count5+1):
    positions.append(i)
limites_final.append(limites[0])
for i in range(1, len(limites)):
    limites_final.append(limites[i]+ limites_final[i-1])

for i in range(1, len(list_positives)):
    print list_positives[i]
    if list_positives[i] != '':
        for x in range(0, len(list_positives[i])):
            list_positives_final.append(int(list_positives[i][x]) + limites_final[i-1])

omegas_final = omegas[0]+omegas[1]
for i in range(2,len(omegas)):
    omegas_final = omegas_final+omegas[i]
seleccion_positiva=0
print lnL_out
for i in range(0, len(lnL_out)):
    if lnL_out[i] > 5:
        seleccion_positiva=1
        break

    
    
if seleccion_positiva != 0 and list_positives_final:
            out= open('/home/devani/PIPE/%s/output.txt' %fold[0], 'w')
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
                out.write('%s' %(omegas_final[i])) 
                out.write('\t')
                for x in range(0, len(list_positives_final)):
                        if list_positives_final[x] == positions[i]:
                            out.write('*')
                            out.write('\t')
                            alarma=1
                            break
                if alarma==0:
                        out.write('-')
                        out.write('\t')
                for x in range(0, len(breakpoints1)):
                    if ((int(breakpoints1[x]))/3) == positions[i]:
                            out.write('*')
                            alarma2=1
                            break
                if alarma2==0:
                        out.write('-')
                out.write('\n')
            out.close()              
                #for i in range(0, len(positions)):
                #    alarma=0
                #    out.write('%s' %(positions[i]))
                #    out.write('\t')
                #    out.write('%s' %(omegas_final[i])) 
                #    out.write('\t')
                #    for x in range(0, len(list_positives_final)):
                #        if list_positives_final[x] == positions[i]:
                #            out.write('*')
                #            alarma=1
                #            break
                #    if alarma==0:
                #        out.write('-')
                #    out.write('\n')
                #out.close()
elif seleccion_positiva != 0 and not list_positives_final:
            out= open('/home/devani/PIPE/%s/output.txt' %fold[0], 'w')            
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
                alarma2=0
                out.write('%s' %(positions[i]))
                out.write('\t')
                out.write('%s' %(omegas_final[i])) 
                out.write('\t')
                for x in range(0, len(breakpoints1)):
                    if ((int(breakpoints1[x]))/3) == positions[i]:
                            out.write('*')
                            alarma2=1
                            break
                if alarma2==0:
                        out.write('-')
                out.write('\n')                
            out.close()
else:
    out= open('/home/devani/PIPE/%s/output.txt' %fold[0], 'w')
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
                alarma2=0
                out.write('%s' %(positions[i]))
                out.write('\t')
                out.write('%s' %(omegas_final[i])) 
                out.write('\t')
                for x in range(0, len(breakpoints1)):
                    if ((int(breakpoints1[x]))/3) == positions[i]:
                            out.write('*')
                            alarma2=1
                            break
                if alarma2==0:
                        out.write('-')
                out.write('\n')
        
out.close()


os.system('rm /home/devani/PIPE/%s/mierda' %fold[0])
os.system('rm /home/devani/PIPE/%s/modelos.txt' %fold[0])
os.system('rm /home/devani/PIPE/%s/otro.txt' %fold[0])


    


#==============================================================================
#         Tenemos dos listas; omegas, con el valor de omega averaged para todos los sitios y positions con la posicion.
#         Van a estar ordenados porque se van introduciendo valores en la lista en el orden de los outs que
#         fue generando el HYPHY
#==============================================================================


    




    
    








