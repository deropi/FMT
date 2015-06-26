#!/usr/bin/env python

protSEQ= open('/home/devani/PIPE/H1.2mammals/protSEQ.fas', 'w')
dnaSEQ= open('/home/devani/PIPE/H1.2mammals/dnaSEQ.fas', 'w')
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
    start = sequence.find('ATG')
    #print start
    #start2=sequence[int(start)+3:].find('ATG')
    #print start2
    #shape= len(sequence[int(start): int(start2)])
    #print shape
    sequencestart = sequence[int(start):]
    start2=sequencestart[3:].find('ATG')
    print len(sequencestart[0:int(start2)])
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


inputsequ = open('/home/devani/PIPE/H1.2mammals/secuencias.fas', 'r')
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


