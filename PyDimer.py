""" 
Takes FASTA sequences and looks for heterodimer interactions. 
 
Creates output .csv with oligo names, sequences, relevant scores, 
and the best alignment for each pair.   
 
Created on Thu Sep 10 19:31:50 2020 
 
@author: ATyler 
"""

# import required modules
import tkinter
from tkinter import filedialog
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import primer3
import pandas as pd

# import fasta sequences to compare
# have user select .fasta file
root = tkinter.Tk()
root.withdraw()  # keeps root window from appearing since no GUI
input_file = filedialog.askopenfilename(
    parent=root,
    title="Select FASTA file to import",
    filetypes=[('FASTA files', '*.fasta')])

# read fasta file into the program
# generate two sequence lists, the original and a list to remove each from
# after completing all pairs with one sequence
static_dict = {rec.id: rec.seq for rec in SeqIO.parse(input_file, "fasta")}
shrink_dict = {rec.id: rec.seq for rec in SeqIO.parse(input_file, "fasta")}
# generate a table to hold each output type
Oligo1 = []
Oligo2 = []
Seq1 = []
Seq2 = []
Tm = []
dG = []
align = []
# sequentially run each sequence in list 1 against all sequences in list 2,
# generating heterodimer data
# iterate over every sequence in the static dict, comparing it to all sequences
# in the "shrink" dict.
for o1 in static_dict.keys():
    for o2 in shrink_dict.keys():
        # ensure the same primer isn't ran against itself
        if static_dict[o1]._data != shrink_dict[o2]._data:
            # calc heterodimer information
            result = primer3.calcHeterodimer(static_dict[o1]._data,
                                             shrink_dict[o2]._data)
            # if structure is found, append relevant values to dict
            if result.structure_found:
                Oligo1.append(o1)
                Oligo2.append(o2)
                Tm.append(result.tm)
                dG.append(result.dg)
                Seq1.append(static_dict[o1]._data)
                Seq2.append(shrink_dict[o2]._data)
                # generate pairwise alignment with no gaps
                aln = pairwise2.align.localxs(static_dict[o1]._data,
                                              shrink_dict[o2].reverse_complement()._data,
                                              -100, -100)
                # collect the alignment information and format as string
                top_aln = aln[0]
                fin_aln = (format_alignment(*top_aln, full_sequences=True))
                # change bottom string back to original sequence
                # convert alignment string to list for iteration
                fin_aln = list(fin_aln)
                # provide key for conversion
                char = {"G": "C", "C": "G", "A": "T", "T": "A"}
                # iterate over string, converting complements after the "alignment"
                # line (which always contains at least one "|")
                bot = False
                for index, item in enumerate(fin_aln):
                    for key, value in char.items():
                        if item == key and bot:
                            fin_aln[index] = value
                        if item == "|":
                            bot = True
                # join corrected string back together
                fin_aln = "".join(fin_aln)
                # append to the relevant dict for storage
                align.append(fin_aln)
    # when each sequence has been tested against all remaining sequences,
    # remove from list 2 so that pairs aren't duplicated
    shrink_dict.pop(o1)
# store each set of sequences and their heterodimer data in the matrix

data = {'Oligo1': Oligo1, 'Oligo2': Oligo2, 'Seq1': Seq1, 'Seq2': Seq2, 'Tm': Tm, 'dG': dG, 'Alignment': align}
df = pd.DataFrame(data)

# export as CSV
export_file_path = filedialog.asksaveasfilename(defaultextension='.csv',
                                                title="Choose filename for results")
df.to_csv(export_file_path, index=False, header=True)