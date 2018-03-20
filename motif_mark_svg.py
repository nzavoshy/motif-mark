#!/usr/bin/python
import re
import argparse
import cairo


parser = argparse.ArgumentParser(description="A program to visualize motifs in a given sequence. The output is a .png file of the location of the motif along the transcript. Showing both introns and exons.")
parser.add_argument("-f", "--file", help='UCSC formatted input sequences, where exons are in capitols and introns are lowercase',required=True,type= str)
parser.add_argument("-m", "--motif", help='a file with the specific motifs to search for within the sequence file,one per line',required=True,type=str) 
args = parser.parse_args()

file=args.file
motif=args.motif

IUPAC_dict={}
IUPAC_dict['R']='[AaGg]'
IUPAC_dict['Y']='[CcTt]'
IUPAC_dict['S']='[CcGg]'
IUPAC_dict['W']='[AaTt]'
IUPAC_dict['U']='[Tt]'
IUPAC_dict['K']='[GgTt]'
IUPAC_dict['M']='[AaCc]'
IUPAC_dict['B']='[CcGgTt]'
IUPAC_dict['D']='[AaGgTt]'
IUPAC_dict['H']='[AaCcTt]'
IUPAC_dict['V']='[AaCcGg]'
IUPAC_dict['N']='[AaCcTtGg]'


def ambig_chr(motif):
'''This function handles ambiguous characters''' #we need a function to turn the iupac names into the actual gene sequences
    save=''
    save_chr=''
    
    for chr in motif:
        chr=chr.upper()
        if chr in IUPAC_dict:
            save+=IUPAC_dict[chr]
            
        else:
            save+='['+chr+chr.lower()+']'
            
            
    return save
with open(motif, 'r') as motif:
    motif_list=[]

    motifs_key=[]
    for motif in motif:

        motif=motif.strip('\n')
        motifs_key.append(motif)
        test=ambig_chr(motif)
        motif_list.append(test)

def make_motifs(cols,n):
'''Function to draw the motifs on the chromosomes we are interested in '''
    context.set_line_width(3)   #this is the location of the motif. #####NEED TO COLOR THE MOTIFS SEPERATELY
    context.set_source_rgb(cols[0],cols[1],cols[2])
    context.move_to(50+loc,135+start_val_mot)        #(x,y)
    context.line_to(50+loc,165+start_val_mot)
    context.stroke()
    context.set_source_rgb(cols[0],cols[1],cols[2])
    context.move_to(1350,80+n)        #(x,y)
    context.line_to(1360,80+n)
    context.stroke()
    #context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, 
    #cairo.FONT_WEIGHT_NORMAL)
    
def make_key(motifs_key):
''' Generates a key for our figure so that we know which motifs have which colors'''
    context.rectangle(1200,50,200,100)
    context.stroke()
    context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, 
    cairo.FONT_WEIGHT_NORMAL)
    context.set_source_rgb(0, 0, 0)
    context.set_font_size(15)
    context.move_to(1210, 70)
    context.show_text('Motif Key')
    n=0
    for key in motifs_key:
    
        context.set_font_size(13)
        context.move_to(1210, 85+n)
        context.show_text(key) 
        n+=15
        context.set_source_rgb(0, 0, 0)

            
#     return(save)
def find_motif(motif_list,line):
''' find the motifs in the sequence file provided'''
    save=[]
    for value in motif_list:
        #print(value)
        for match in re.finditer(value,line): #find the motif in the line
            #print(match)
            s=match.start() #find the starting position of the match
            #e=match.end()
            save.append((value,s))
            #save.append(s)
            #save.append(e)
            
    return(save)
def positions(line,motif_list): # we need to find the starting and ending position in order to actually map it
'''Extracts important info from the header. Including: chromsome, start and end position, and name'''
#part of this can be taken out of the header of the FASTA 
    if line.startswith('>'):
        chrome=line.split(' ')[0][1:]
        header= line.split(' ')[1]
        match=re.search(r'([\w]+)\:([0-9]+)\-([0-9]+)',header) #grab the header and separate it into groups
        chrom_num=match.group(1) #chromosome
        start_pos=match.group(2) #starting position
        end_pos=match.group(3) #ending position
        return(chrome,chrom_num,int(end_pos)-int(start_pos)) #return a toupple of this information

        
        #return(s,e)
def find_exon(file,motif_list):

    '''This function will convert the file, minus the header line, into a string and then search for the exons location'''##this will show you where the exon is in your file with the beginning and ending coordinates. Need to minus 
#one because it will actually tell you where the first lowercase begins 
    sequence='' #empty sequence
    options=[] #list of the things to return
    first_line=True #this is necessary
    for line in file: #go through each line in the file
        if first_line and line.startswith('>'):#if its the first line, extract header
        #if line.startswith('>'):
            header=positions(line,motif) #positions function grabs the header line
            first_line=False
        elif not first_line and line.startswith('>'):
            
            #print(line)
            save_me=find_motif(motif_list,sequence)
            match=re.search(r'[ACTG]+',sequence)
            #print(sequence)
            #print(match)
            s=match.start()
            e=match.end()
            save=(s,e,header)
            
            #print(save_me)
            options.append(save)
            options.append(save_me)
            header=positions(line,motif)
            sequence=''
        else:
            line=line.strip('\n')
            #print(line)
            sequence+=line
            #print(sequence)
    #print(sequence)
    match=re.search(r'[ACTG]+',sequence)
    s=match.start()
    e=match.end()
    save=(s,e,header)
    options.append(save)
    save_me=find_motif(motif_list,sequence)
    options.append(save_me)
    #print(save_me)
    return(options)
with open(file, 'r') as file:
    exon=find_exon(file,motif_list)

            
    
    
    #print(exon)
start_val=100
start_val_mot=100
n=0
m=0

import random
cols={}
nums=[]
for x in motif_list: #make a dictionary of colors to use for however many motifs you have. 
    nums=[]
    for y in range(3):
        num1=random.random() 
        nums.append(num1)

    cols[x]=nums
#print(cols)
#print(cols)
width, height = 1500, 1000 #initiate the drawing
#create the coordinates to display your graphic, desginate output
surface = cairo.SVGSurface("motifs.svg",width, height) 
#create the coordinates you will be drawing on (like a transparency) - you can create a transformation matrix
context = cairo.Context(surface)
for value in exon:
    #print(value)

    if type(value)==tuple:
        header=value
        title=header[2][0]
        chrome=header[2][1]
        trans_len=header[2][2]
        exon_start=header[0]
        exon_end=header[1]
        print(trans_len)
        print(exon_start)
        print(exon_end)
    
        context.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, 
        cairo.FONT_WEIGHT_NORMAL)
        context.set_source_rgb(0, 0, 0)
        context.set_font_size(13)
        context.move_to(20, 30+start_val)
        context.show_text(str(title))
        context.move_to(20, 45+start_val)
        context.show_text(str(chrome))
        
        context.set_line_width(3)    #this makes the straight line(which will be the gene)
        context.set_source_rgb(0, 0, 0)
        context.move_to(50,150+start_val)        #(x,y)
        context.line_to(50+trans_len,150+start_val)
        context.stroke()


        context.set_source_rgb(0, 0, 0)
        context.rectangle(exon_start,140+start_val,exon_end-exon_start,20)
        context.fill()

        n+=1
        start_val+=200
    else:
        #start_val_mot=100
        first_motif=motif_list[0]

        n=0
        for item in value:
            motif=item[0]
            loc=item[1]
            if motif==first_motif:
                #print(motif)
                x=cols[motif]
                make_motifs(x,n)

            else:
                n+=15
                #print(motif+'not first motif')
                x=cols[motif]
                make_motifs(x,n)
                first_motif=motif

        context.set_source_rgb(0, 0, 0)
            
            
        start_val_mot+=200

make_key(motifs_key)

surface.finish()
        


