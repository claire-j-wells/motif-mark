#!/usr/bin/env python

import argparse
import cairo
import math

def get_args():
    parser = argparse.ArgumentParser(description="A script that finds motifs for up to 10 sequences and up to 5 motifs")
    parser.add_argument("-f", help="input fa file", required= True)
    parser.add_argument("-m", help="input motifs_file", required = True)
    
    return parser.parse_args()

args = get_args()
file = args.f 
motif = args.m

#FUNCTIONS

def oneline_fasta(file_input,first_output):
    first_line=True
    with open(first_output,"w") as op: ##to get the files on a single line
        with open(file_input,"r") as file:
            for line in file:
                line=line.strip('\n')
                if line.startswith(">"):
                    if first_line==True:
                        op.writelines([line,'\n'])
                        first_line=False
                    else:
                        op.writelines(['\n',line,'\n'])
                else:
                    op.writelines(line)

#We want to split the exons and introns for the drawing to know whether we want a line or a box. 
def exon_intron_split(sequence):
    uplow = sequence[0].isupper() # boolean to check whether the first letter is uppercase or not
    splitted = []
    new_word = ""
    position_counter = 0
    for c in sequence:
        if c.isupper() == uplow:
            new_word += c
            position_counter+=1
        else:
            if c.isupper()!=True:
                splitted.append(("E",new_word,len(new_word),position_counter-len(new_word)))
            else:
                splitted.append(("I",new_word,len(new_word),position_counter-len(new_word)))
            new_word = c
            position_counter+=1
            uplow = not uplow
    if c.isupper()==True:
                splitted.append(("E",new_word,len(new_word),position_counter-len(new_word)))
    else:
                splitted.append(("I",new_word,len(new_word),position_counter-len(new_word)))
      # bc no change at end of string
    return(splitted)

def motif_permutations(motif):
    wobbles = {
    "Y":["T","C"],
    "W":["A","T"],
    "S":["C","G"],
    "M":["A","C"],
    "K":["G","T"],
    "R":["A","G"],
    "B":["C","G","T"],
    "D":["A","G","T"],
    "H":["A","C","T"],
    "V":["A","C","G"],
    "N":["A","T","C","G"]
}
    motif = motif.upper()
    all_poss = [[]]
    #for each char 
    for character in motif:
        #if normal sequence add to all our curr poss
        if character not in wobbles:
            all_poss=[item+[character] for item in all_poss]
            
        else:
            #get out possible wobbles
            poss=wobbles[character]
            #count our possible wobbles
            num_poss=len(poss)
            #determine how many permutations we had before this char
            one_poss=len(all_poss)
            #track which wobble we should use
            wobble_index=0
            #determine how many perms this wobble will produce
            all_poss=all_poss*num_poss
            #go through all the permuations
            for i in range(len(all_poss)):
                #if we have exhausted the possibilites for this wobble
                if i>=one_poss:
                    #move to the next one
                    one_poss+=one_poss
                    wobble_index+=1
                #add our current wobble to our permuations
                all_poss[i]=all_poss[i]+[poss[wobble_index]]
    #Makes a string
    all_poss=list(set(["".join(item) for item in all_poss]))
    return(all_poss)


#CLASSES

class Gene:
    name = None
    sequence = None           
    sequence_len = None
    split_exons_introns = None
    sequence_single_case = None
    #Gene object must take in a header and a sequence 
    def __init__(self, header, the_sequence):
        #Grab the name of the gene 
        self.name = header.strip("\n")
        self.sequence = the_sequence
#For the introns and exons              
        self.sequence_len = len(the_sequence)
        #self.intron_exon = re.split('([A-Z]+)',sequence) #splits the sequence by exons and introns
        self.split_exons_introns = exon_intron_split(the_sequence)
#For the motifs
        self.sequence_single_case = the_sequence.upper() 

class Motifs:
    gene_name = None
    name = None 
    length = None            
    motif_var = None
    position = None
    def __init__(self,motif): 
        self.gene_name = []
        self.name = motif.strip("\n")              
        self.length = len(self.name)            
        self.motif_var = motif_permutations(motif)
        self.position = []

class Storage_Bin:
    gene_name = None
    motif = None
    position = None 
    color = None
    mot_length = None
    motif_orig = None
    def __init__(self, gene_name, motif, position, color, mot_length, motif_orig):
        self.gene_name = gene_name
        self.motif = motif
        self.position = position
        self.color = color
        self.mot_length = mot_length
        self.motif_orig = motif_orig


#ESTABLISH A LIST OF GENE OBJECTS

#for future use, we want to split the file input into variables to fstring her later
file_name = file
#split by period to only get the prefix
split = file_name.split(".")
#index out the prefix and assign it to a variable
file_pre = split[0]

#run one_line fasta 
oneline_fasta(file, f'{file_pre}_2.fa')

#pass one line fasta output into for look to make gene objects
with open(f'{file_pre}_2.fa', "r") as the_file:

    #store gene objects in a list to iterate through later
    gene_obj_list = []
    for line in the_file:
        if line.startswith(">"):
            head = line
        else:
            seq = line
            g=Gene(head,seq)
            gene_obj_list.append(g)


#ESTABLISH MOTIF OBJECTS

with open(motif, "r") as mot:

    #motif_list is a list of motif OBJECTS to iterate through later
    motif_list = []
    for line in mot:
        line = line.strip("\n")
        line = line.replace("U", "T")
        m=Motifs(line)
        motif_list.append(m)

#list of colors to index through during the drawing process 
colors_list= [(51, 180, 10),(57, 48, 73),(218, 21, 177),(80, 213, 255),(244, 178, 65)]

#ESTABLISH STORAGE BIN OBJECTS FOR DRAWING

#There should be a storage bin object for every single motif that is found in the FASTA file 
def grab_position_motifs(gene_obj_list,motif_list):
    '''This function takes in two lists of objects (gene objects and motif objects)
    and it passes the information from these object into a new storage bin object for drawing '''
    storage_bin_list = []
    for gene_object in gene_obj_list:
        gene_name = gene_object.name
        for j,letter in enumerate(gene_object.sequence_single_case):
            #print(letter)
            for i, objects in enumerate(motif_list):
                motif_orig = objects.name
                motif_name = objects.motif_var
                mot_length = objects.length
                motif_checker = gene_object.sequence_single_case[j:j+objects.length]
                if motif_checker in objects.motif_var:
                    position = j
                    color = colors_list[i]
                    match = Storage_Bin(gene_name, motif_name, position, color, mot_length, motif_orig)
                    storage_bin_list.append(match)
                    #print(sub_list) 
        #print(motif_loc)
    return(storage_bin_list)

#assign the outputted list from the function above to a variable. 
storage_bin_list=grab_position_motifs(gene_obj_list,motif_list)

gene_lengths_list = []
for gene_objects in gene_obj_list:
    gene_lengths_list.append(gene_objects.sequence_len)
cairo_width = max(gene_lengths_list)
cairo_height = len(gene_lengths_list)*300 + 500
WIDTH,HEIGHT = cairo_width, cairo_height
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, HEIGHT)
context = cairo.Context(surface)

x_start = 0
y_start = 160
y_increase = 2


context.rectangle(0, 0, cairo_width, cairo_height)
context.set_source_rgb(1, 1, 1)
context.fill()
context.set_font_size(25)
context.select_font_face(
        "Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)


for gene_object in gene_obj_list:
    context.set_source_rgb(0, 0, 0)
    context.move_to(x_start,y_start)# y_increase)
    context.line_to(gene_object.sequence_len,y_start)
    context.stroke()
    context.move_to(20, y_start-50)
    context.text_path(gene_object.name)
    for tuples in gene_object.split_exons_introns:
        #print(tuples[1])
        if tuples[0][0] == "E":
            x_rect = tuples[3]
            y_rect = y_start-20
            w_rect = tuples[2]  
            h_rect = 40
            context.rectangle(x_rect, y_rect, w_rect, h_rect)
            context.fill()
        else:
            continue
    for box in storage_bin_list:
        #print(box.motif_orig)
        if box.gene_name == gene_object.name:
            #print(box.motif_orig)
            #print(box.color)
            context.move_to(x_start,y_start)
            context.set_source_rgba(box.color[0]/255,box.color[1]/255, box.color[2]/255, 0.6)
            x_rect = box.position
            y_rect = y_start-40
            w_rect = box.mot_length  
            h_rect = 80
            context.rectangle(x_rect, y_rect, w_rect, h_rect)
            context.fill()       
    y_start=y_start+300

#BUILD THE LEGEND

#this allows the legend to move around and start drawing relative 
#to how many genes are present in the FASTA file

legend_start_y = y_start - 100
move_text_y = y_start - 80

#Build a dictionary with each motif and its respective color 
motif_color_dict = {}
for box_object in storage_bin_list:
    motif_color_dict[box_object.motif_orig] = box_object.color

#Seperately add an Exon entry because Exons are not stored as an object in storage_bin
motif_color_dict.update({"Exon": (0,0,0)})

#Due to issues with opacity in the lends, Exons are drawn first outside the for-loop
for key in motif_color_dict:
    if key == "Exon":
        context.set_source_rgb(motif_color_dict[key][0]/255,motif_color_dict[key][1]/255, motif_color_dict[key][2]/255)
        context.rectangle(legend_start_x, legend_start_y, 70, 20)
        context.fill()
        context.move_to(move_text_x,move_text_y)
        context.set_source_rgb(0,0,0)
        context.text_path(key)
        context.stroke()
        legend_start_y = legend_start_y + 50
        move_text_y = move_text_y + 50
#If the key does not equal Exon, then draw the legend normally AND include opacity 
    else:       
        legend_start_x = 50
        move_text_x = 150
        context.set_source_rgba(motif_color_dict[key][0]/255,motif_color_dict[key][1]/255, motif_color_dict[key][2]/255, 0.6)
        context.rectangle(legend_start_x, legend_start_y, 70, 20)
        context.fill()
        context.move_to(move_text_x,move_text_y)
        context.set_source_rgb(0,0,0)
        context.text_path(key)
        context.stroke()
        legend_start_y = legend_start_y + 50
        move_text_y = move_text_y + 50

#Write to surface with the correct prefix name 
surface.write_to_png(f'{file_pre}.png')