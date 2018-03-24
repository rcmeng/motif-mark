#!/usr/bin/env python

####Argeprase Setup####


#Importing argparse so that the program can be run from the command line, cairo so that the location of the motifs with respect to the introns
#and exons can be visualized, and finally random for random number generation for the colors of our output.
import argparse
import cairo
import random


#Desiginating the required arguments for this program. These include the input fasta file (ideally from the UCSC browser
#so that intron and exon sequences are defined) and a text file of motifs. 

parser = argparse.ArgumentParser(description="Given a fasta file from the UCSC genome browser (desiginating introns and exons)and a file of motifs (using IUPAC designation), this program produces a SVG file visualizing the location of these motifs with respect to the introns and exon of the sequence.")
parser.add_argument("-f", "--sequence_path", help="Requires the path to the fasta file that needs to be parsed and visualized.")
parser.add_argument("-m", "--motif_path" ,help="Requires the path to the motif text file that will be used to visualize the motifs with respect to the introns and exons")
args = parser.parse_args()
fasta_file= args.f
motif_file = args.m


####Python Function Definitions####


#Defining a dictionary for later use in the program with the IUPAC base designations
IUPAC_dict = {"A": "A", "C":"C", "G":"G", "T":"TU", "U":"UT",
              "R": "AG", "Y":"TCU", "S":"CG", "W":"ATU",
              "K":"GTU", "M":"AC", "B":"CGTU", "D":"AGTU",
              "H":"ACTU", "V":"ACG", "N":"ATCGU"}
			  
			  
def sequence_parser(sequence):
'''Takes an entry from a FASTA file and parses it into two introns and an exon where applicable'''
	#Setting empty strings for the two introns and the exon from the fasta file
	intron1 = ""
	exon = ""
	intron2 = ""
	exon_presence = False
	#Looping through each fasta entry distinguishing between the two introns and exons, setting each empty string to the actual introns and exons
	for character in line:
		if character.isupper():
			exon += character
			exon_presence = True
		elif exon_presence = True:
			intron2 += character
		else:
			intron1 += character
	return [intron1, exon, intron2]
	
def fasta_parser(fasta_file):
'''Takes an input fasta file and seperates each entry into the header and sequence lines and then sets those as key-value pairs in a dictionary'''
    #Looping through the fasta file, setting an empty dictionary and empty strings for the header and sequence. Also
	#setting a variable to check to see if the entry is complete to store the value and move onto the next header/sequence pair.
	with open(fasta_file)as fh:
        complete_entry = False
	    fasta_dict = {}
        header = ""
        sequence = ""
		#Looping through each line setting each header as the key in our empty dictionary, and each sequence line as the associated value
		for line in fh:
            line = line.strip()
			if line.startswith('>'):
                if complete_entry == True:
                    fasta_dict[header] = sequence_parser(sequence)
                header = line
                sequence = ""
            else:
                sequence += line
                complete_entry = True
	fasta_dict[header] = sequence_parser(sequence)
    return fasta_dict


def motif_parser(motif_file):
'''Iterates through a file of motifs and creates a list of possible characters for each motif.'''
    #Setting variables for an empty motif list and empty motif dictionary
	motif_list = []
    motif_dict = {}
	#Iterating thourgh the motif file stripping the new line character, making all the characters uppercase, and appending them to our 
	#aforementioned motif list
    with open(motif_file)as fh:
        for line in fh:
            line = line.strip()
            line = line.upper()
            motif_list.append(line)
	#Looping through the motif list and setting each motif as a key and the variety of bases that it can represent to the corresponding value.
    for motif in motif_list:
        chars = []
        for characters in motif:
            chars.append(base_dict[characters])
        motif_dict[motif] = chars
	return motif_dict

def get_positions(header, sequence, motif_dict):
'''Function that creates a dictionary to store the position of each of the motifs passed to it'''
	#Setting an empty dictionary for position
	positions = {}
	#Iterating through our aformentioned motif dictionary setting chars as the value of the motif and a variable for the length of the motif.
    for motif in motif_dict.keys():
        chars = motif_dict[motif]
        motiflength = len(motif)
		#Nested for loop that creates the position dictionary which has the header plus motif as the key, and the location of the motifs
		#as the value pair.
        for i in range(0,len(sequence) - motiflength):
            pattern = sequence[i:i + motiflength]
			count = 0
            motif_match = True
            for c in pattern:
				if c not in chars[count]:
                    motif_match = False
                count += 1
            if match == True:
				if header + "_" + motif in positions:
                    positions[header + "_" + motif].append(i)
                else:
                    positions[header + "_" + motif] = [i]
	return positions


def motif_finder(fasta_file, motif_file):
'''Finds positions of all motifs and returns them in a list format'''
    #Utilizing our previous functions to parse the fasta and motif files, creating new fasta and motif dictionaries in 
	#addition to initalizing an empty list of fasta positions.
	fasta_dict = fasta_parser(fasta_file)
    motif_dict = motif_parser(motif_file)
    fasta_position_list = []
	#Iterating through all of the items in the fasta dicttionary to make a list of the positions of the motifs in all of the introns and exons.
    for item in fasta_dict.items():
        seq = "".join(fasta_dict[item[0]])
		pos = get_positions(item[0], seq, motif_dict)
		pos["intron1"] = len(fasta_dict[item[0]][0])
        pos["exon"] = len(fasta_dict[item[0]][1])
        pos["intron2"] = len(fasta_dict[item[0]][2])
		fasta_pos_list.append(pos)
	return fasta_pos_list

####PyCairo drawing functions####

def draw_exon(intron1 ,exon_length, intron2, context, start):
'''Draws the lines for each of the two introns and exons in our drawing'''
    #Setting a line width of 1 and source colors
	context.set_line_width(1)
    context.set_source_rgb(0,0,0)
	#Moving to the starting position and making a stroke for the line itself
    context.move_to(start[0], start[1])
    context.line_to(start[0] + intron1 + exon_length + intron2, start[1])
    context.stroke()
	#Drawing the exon itself by utilizing a larger line width
	context.set_line_width(15)
    context.move_to(start[0] + intron1, start[1])
    context.line_to(start[0] + intron1 + exon_length, start[1])
    context.stroke()

def draw_motifs(pos_list, context, start, motif):
'''Draws the motifs onto the two introns exon. Also returns the color of each individual motif'''
    #Randomly generating colors
	red = random.random()
    green = random.random()
    blue = random.random()
	#Setting the line width and source colors for our motifs
    context.set_line_width(10)
    context.set_source_rgb(red, green, blue)

	#Iterating through our position list and drawing the motifs in their designated location
    for pos in pos_list:
        context.move_to(start[0] + pos, start[1])
        context.line_to(start[0] + pos + len(motif), start[1])
        context.stroke()

        context.move_to(start[0] + pos, start[1] - 20)
        context.show_text(str(pos))

    return [red, green, blue]

def draw_legend(red, green, blue, start, drop, motif):
'''Draws the legend with the colors of the motifs'''
    context.set_source_rgb(red, green, blue)
    context.move_to(start[0] - 50, start[1] + drop)
    context.show_text(motif)

def draw_title(start, name):
'''Draws the name of each gene next to the legend'''
    context.set_source_rgb(0,0,0)
    context.move_to(start[0] - 100, start[1])
    context.show_text(name)


#Creating a list of dictionaries of motif positions for each fasta entry
all_pos = motif_finder(fasta_file, motif_file)

#Establishig the PyCairo surface and setting the position to start
surface = cairo.SVGSurface("plot.svg", 1850, 1200)
context = cairo.Context(surface)
context.set_line_width(1)
start = [300,200]

#Using our aforementioend draw exon function to draw the introns and exons
for f in all_pos:
    draw_exon(f["intron1"], f["exon"], f["intron2"], context, start)
	#Centering the text of the legend
    drop = -20 
    for m in f.items():
        if m[0][0] == ">":
            motif = m[0].split("_")
            motif = motif[1]
            gene = m[0].split(" ")
            gene = gene[0][1:]
			color = draw_motifs(m[1], context, start, motif)
            draw_title(start, gene)
            draw_legend(color[0], color[1], color[2], start, drop, motif)
            drop = drop + 10
	start = [start[0], start[1] + 100]