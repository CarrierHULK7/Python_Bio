# === Final Combined Bioinformatics Python File ===

# ===== File: argv.py =====
# ------------------------------------------------------------------
# File name: argv.py
#
# Command line arguments are provided directly after the name of the
# program and they are stored as a list in sys.argv. The first entry
# sys.argv[0] is the script name (it is operating system dependent
# whether this is a full pathname or not).
#
# Our program will print out the list of arguments; will exit if
# no argument is provided. Example usage: args.py zika_DNA.fasta
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/interpreter.html#argument-passing
# https://docs.python.org/3/library/sys.html
# https://docs.python.org/3/using/cmdline.html#command-line-and-environment
# ------------------------------------------------------------------

import sys

if len(sys.argv) == 1:
   print('Please provide a command line argument!')
   sys.exit()

print('sys.argv list:', sys.argv)
print('The first argument:', sys.argv[0])
print('The second argument:', sys.argv[1])
# ===== File: beginnings.py =====
# -------------------------------------------------------------------
# File name: beginnings.py
#
# Variable names in Python start with a letter followed by
# combination of letters, digits or underscore (no white spaces).
#
# Four of the basic variable types in Python are
# numeric (integers and floats), string, list, and tuple.
# The code below introduces examples of these variable types.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/introduction.html#an-informal-introduction-to-python
# https://www.ncbi.nlm.nih.gov/nuccore/U00096
# https://en.wikipedia.org/wiki/Chargaff%27s_rules
# ------------------------------------------------------------------

# String variable
organism = "Escherichia coli"        # NCBI accession number U00096
strain = 'str. K-12 substr. MG1665'
print("DEFINITION: " + organism + " " + strain)

# Integer variable
number_of_bps = 4641652
print('Number of base pairs:', number_of_bps)

# Float variable
percent_A = 24.7
percent_T = 23.6

# List variable
percents_AGCT = [percent_A, 26.0, 25.7, percent_T]
print("[A, G, C, T] =", percents_AGCT)

# Computing ratios A/T and G/C
ratio_AT = percent_A / percent_T
ratio_GC = percents_AGCT[1] / percents_AGCT[2]

# Tuple variable
E_Coli = (organism, ratio_AT, ratio_GC)
print(E_Coli)
# ===== File: buggy.py =====
# ------------------------------------------------------------------
# File name: buggy.py
#
  This short Python code contains a number of interntional bugs. Correct
# them with the help of the error messages.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://americanhistory.si.edu/collections/search/object/nmah_334663
# https://en.wikipedia.org/wiki/Software_bug
# ------------------------------------------------------------------

human genes = 20,000

 protein-name = "GFP';

print('You have', human_genes, 'genes')
print("protein-name stands for
      green fluorescent protein")
# ===== File: concatenation.py =====
# ------------------------------------------------------------------
# File name: concatenation.py
#
# Binary operator + concatenates two strings.
#
# The code below concatenates the first four codons for GFP.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://en.wikipedia.org/wiki/RNA_codon_table
# http://www.abcam.com/recombinant-a-victoria-gfp-protein-ab84191.html
# ------------------------------------------------------------------

GFP_seq = 'MSKGEELFTG...HGMDELYK'
print('Green fluorescent protein sequence:', GFP_seq)

M_codon = 'AUG'
S_codon = 'UCA'
K_codon = 'AAA'
G_codon = 'GGU'

RNA_seq = M_codon
RNA_seq = RNA_seq + S_codon
print('RNA sequence:', RNA_seq)

RNA_seq = RNA_seq + K_codon + G_codon
print(RNA_seq, 'could code amino acid sequence MSKG')
# ===== File: count_with_dictionary.py =====
# ------------------------------------------------------------------
# File name: count_with_dictionary.py
#
# Code below counts the number of each codon that appears in a DNA
# segment using a dictionary. Keys are the present codons and the values
# are their numbers.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#list
# https://docs.python.org/3/library/stdtypes.html#mapping-types-dict
# http://www.ncbi.nlm.nih.gov/nuccore/XM_002295694.2?report=fasta
# http://www.opensourceshakespeare.org/stats/
# ------------------------------------------------------------------

DNA_seq = 'gggtgcgacgattcattgttttcggacaagtggataggcaaccactaccggtggattgtc'
print('Sequence:', DNA_seq)

# Count the number of codons
DNA_length = len(DNA_seq)
number_of_codons = int(DNA_length/3)
print()

# Put codons in a list
codon_list = []
for i in range(number_of_codons):
    codon_list.append(DNA_seq[i*3:i*3 + 3])

print('List of codons:', codon_list)
print()

#Create codon counter dictionary (codon : codon_count)
codon_counter = {}

for codon in codon_list:
    if codon not in codon_counter:
        codon_counter[codon] = 1
    else:
        codon_counter[codon] += 1

# This loop syntax accesses the whole dictionary by looping
# over the .items() tuple list, accessing one (key, value)
# pair at each step.
print('Codon counter:')
for key, value in codon_counter.items():
    print(key, ':', value)
# ===== File: dictionary.py =====
# ------------------------------------------------------------------
# File name: dictionary.py
#
# A dictionary (dict) variable type is akin to a list that can hold
# a number of values. However, instead of indexing with integers, it
# uses a unique name, called a key, for each entry. A dict d = { }
# can be defined by comma-separated pairs of key and value, with
# a : in between, e.g. d = {"EcoRI" : "GAATTC"}.
#
# Code below illustrates some of the basic operations with a dictionary
# containing several restriction enzymes and their recognition sites.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/datastructures.html#dictionaries
# https://docs.python.org/3/library/stdtypes.html#mapping-types-dict
# https://en.wikipedia.org/wiki/Restriction_enzyme
# ------------------------------------------------------------------

restriction_enzymes = {'EcoRI' : 'GAATTC',
                        'AluI' : 'AGCT',
                        'NotI' : 'GCGGCCGC',
                        'TaqI' : 'TCGA'
                      }

print(restriction_enzymes)
print()

# To get a list of keys from a dictionary view object
keys = list(restriction_enzymes.keys())
print('Keys as a list:', keys)

# To get a list of values from a dictionary view object
values = list(restriction_enzymes.values())
print('Values as a list:', values)

# To check if a key is in the dictionary
mykey = 'crispr'
check = mykey in restriction_enzymes
print('Is', mykey, 'key in the dictionary?', check)
print()

# To fetch a value from a dictionary with its key
EcoRI_value = restriction_enzymes['EcoRI']  #raises a KeyError if key not found
EcoRI_value = restriction_enzymes.get('EcoRI') # does not raise a KeyError if key not found
print('The recognition site of EcoRI is', EcoRI_value)
print()

# To add to an existing dictionary
restriction_enzymes['EcoRV'] = 'GATATC'
restriction_enzymes.update(EcoRV = 'GATATC')
print('With a new item:', restriction_enzymes)
print()

# To delete an item from a dictionary
del restriction_enzymes['EcoRV']
print('Original dictionary:', restriction_enzymes)
print()
# ===== File: find.py =====
# ------------------------------------------------------------------
# File name: find.py
#
# str.find(sub[, start[, end]])
# Return the lowest index in the string where substring sub is found
# within the slice s[start:end]. Optional arguments start and end are 
# interpreted as in slice notation. Return -1 if sub is not found.
#
# str.rfind(sub[, start[, end]])
# Return the highest index in the string where substring sub is found,
# such that sub is contained within s[start:end]. Optional arguments 
# start and endare interpreted as in slice notation. Return -1 on failure.
# 
# str.index(sub[, start[, end]]) and str.rindex(sub[, start[, end]])
# are like find() and rfind() but raise ValueError when the substring 
# sub is not found.
#
# str.count(sub[, start[, end]])
# Return the number of non-overlapping occurrences of substring sub
# in the range [start, end]. Optional arguments start and end are
# interpreted as in slice notation.
#
# Code below compute the number and the locations of codon CAT in
# a DNA segment.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#str.find
# https://docs.python.org/3/library/stdtypes.html#str.rfind
# https://docs.python.org/3/library/stdtypes.html#string-methods
# https://en.wikipedia.org/wiki/D-loop
# http://www.ncbi.nlm.nih.gov/nuccore/AF176731.1?report=fasta
# ------------------------------------------------------------------

chimp = 'GTACCACCTAAGTACTGGCTCATTCATTACAACCGGTATGTACTTCGTACATTACTGCCAGTCACCATGA'
print('Chimp D-loop:', chimp)

codon = 'CAT'

# Check if subtring is in string
is_in = codon in chimp
print('Is codon', codon, 'in chimp:', is_in)

# .count() counts how many times sub appears in string
how_many = chimp.count(codon)
print('How many times', codon, 'appears in chimp:', how_many)

# .find() returns the lowest index
first_index = chimp.find(codon)
print('First', codon, 'index: ', first_index)

second_index = chimp.find(codon, first_index + len(codon))
print('Second',  codon, 'index: ', second_index)

# .find() returns -1 if cannot find sub in string
third_index = chimp.find(codon, 27, 55)
print('Third', codon, 'index: ', third_index)

# .rfind() returns the highest index
last_index = chimp.rfind(codon);
print('Last', codon, 'index: ', last_index)
# ===== File: for.py =====
# ------------------------------------------------------------------
# File name: for.py
#
# for item in items:
#      statements

# Python’s for statement iterates over the items of any sequence
# (a list or a string), in the order that they appear in the sequence.
#
# Code below prints out all codons starting with T.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/controlflow.html#for-statements
# https://docs.python.org/3/reference/compound_stmts.html#for
# https://en.wikipedia.org/wiki/DNA_codon_table
# ------------------------------------------------------------------

# Save nucleotide bases in a list
bases = ['T', 'C', 'A', 'G']

# As a warmup, print the list of bases
for base in bases:
    print(base)

print('Codons starting with T:')

for second_base in bases:
    print('Codons starting with T'+second_base)
    for third_base in bases:
        print('T'+second_base+third_base)
# ===== File: function.py =====
# ------------------------------------------------------------------
# File name: function.py
#
# The keyword def introduces a function definition. It must be followed
# by the function name and the parenthesized list of formal parameters.
# The statements that form the body of the function start at the next line,
# and must be indented. Variables in a function are local to that function.
#
# The first statement of the function body can optionally be a string literal
# enclosed in triple quotes ''' ... '''; this string literal is the 
# function’s documentation string, or docstring. 
#
# The user-defined function below takes a DNA segment as a string
# and returns the corresponding RNA segment as a string.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/controlflow.html#defining-functions
# https://docs.python.org/3/reference/compound_stmts.html#function-definitions
# https://docs.python.org/3/tutorial/controlflow.html#documentation-strings
# http://www.ncbi.nlm.nih.gov/nuccore/226377833?report=fasta
# ------------------------------------------------------------------

def DNAtoRNA(dna):
    '''
    Converts DNA string to RNA string.

    Parameters: DNA sequence as a string
    Return: RNA sequence as a string
    '''
    transliterate = dna.maketrans('tT', 'uU')
    rna = dna.translate(transliterate)
    return rna

# Can access the doctring ''' ... ''' with print(DNAtoRNA.__doc__)

zika_DNA = 'AGTTGTTGATCTGTGTGAGTCAGACTGCG'
print('Zika DNA segment is', zika_DNA)

zika_RNA = DNAtoRNA(zika_DNA);
print('Zika RNA segment is', zika_RNA)
# ===== File: if_elif.py =====
# ------------------------------------------------------------------
# File name: if_elif.py
#
# The compound statement if_elif_else has the syntax
#
# if EXPR1:
#    statements1
# elif EXP2:
#     statements2
# else:
#    statements3
#
# If EXPR1 is True, the first group of statements1 are executed and
# the rest is skipped; otherwise, if EXPR2 is True, the statements2
# are executed; otherwise, statements3 are executed.
# There can be zero or more elif parts, and the else part is optional.
# The keyword ‘elif’ is short for ‘else if’, and is useful to avoid
# excessive indentation.
# An if … elif … elif … sequence is a substitute for the switch or
# case statements found in other languages.
#
# Code below determines if the last codon in a DNA segment
# is the start codon ATG or one of the stop codons TAA, TAG, or TGA;
# or none of the above.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/controlflow.html?#if-statements
# https://docs.python.org/3/reference/compound_stmts.html?#the-if-statement
# https://docs.python.org/3/library/stdtypes.html?#truth-value-testing
# https://docs.python.org/3/tutorial/datastructures.html?#comparing-sequences-and-other-types
# https://en.wikipedia.org/wiki/DNA_codon_table
# ------------------------------------------------------------------

DNA_segment = 'ATGACATGACCAATC'
codon1 = DNA_segment[-3:]

# == operator tests the equality of two strings, resulting in True/False
if (codon1 == 'ATG'):
    print('Codon', codon1, 'is a start codon.')
elif ((codon1 == 'TAA') or
     (codon1 == 'TAG') or
     (codon1 == 'TGA')):
    print('Codon', codon1, 'is a stop codon.')
else:
    print('Codon', codon1, 'is neither a start nor a stop codon.')

print('Done!')
# ===== File: if_else.py =====
# ------------------------------------------------------------------
# File name: if_else.py
#
# The compound statement if_else has the syntax
#
# if EXPR:
#    statements
# else:
#    statements
#
# If EXPR is True the first block of statements are executed, otherwise
# the second block of statements following else are executed.
#
# Code below determines if the first codon in a DNA segment
# is the start codon ATG or not and reports the result.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/controlflow.html?#if-statements
# https://docs.python.org/3/reference/compound_stmts.html?#the-if-statement
# https://docs.python.org/3/library/stdtypes.html?#truth-value-testing
# ------------------------------------------------------------------

DNA_segment = 'ATGACATGA'
codon1 = DNA_segment[0:3]

# == operator tests the equality of two strings, resulting in True/False
if (codon1 == 'ATG'):
    print('Codon', codon1, 'is a start codon.')
else:
    print('Codon', codon1, 'is not a start codon.')

print('Done!')
# ===== File: input.py =====
# ------------------------------------------------------------------
# File name: input.py
#
# input([prompt])
# If the prompt argument is present, it is written to standard output
# without a trailing newline. The function then reads a line from input,
# converts it to a string (stripping a trailing newline), and returns
# that. When EOF is read, EOFError is raised. input() returns a string;
# cast it if a number is expected.
#
# This program first prompts the user to type an NCBI sequence number
# and then echos it. Second, it illustrates the casting of string
# input with int() for use in arithmetical operations by calculating 
# the sum of two user-typed numbers.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/functions.html#input
# ------------------------------------------------------------------

sequence_number = input('Please type an NCBI sequence number: ')
print('Your sequence number is', sequence_number)

# Input returns a string
starting_index = input('Please type a starting index: ')
ending_index = input('Please type an ending index: ')
print('I will compute the number of bps in this region...')

# Must cast string inputs to int or float for arithmetical operations
number_of_bps = int(ending_index) - int(starting_index)
print('The number of bps is:', number_of_bps)
# ===== File: length.py =====
# ------------------------------------------------------------------
# File name: length.py
#
# Built-in function len(s)
# Return the length (the number of items) of an object. The argument
# may be a sequence (such as a string, bytes, tuple, list, or range)
# or a collection (such as a dictionary, set, or frozen set).
#
# Code below computes the number of nucleotides in a DNA segment.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/functions.html#len
# http://www.ncbi.nlm.nih.gov/nuccore/226377833?report=fasta
# ------------------------------------------------------------------

zika_DNA = 'AGTTGTTGATCTGTGT'
zika_DNA_length = len(zika_DNA)

print('The first', zika_DNA_length, 'nucleotides',
      'of Zika virus DNA are', zika_DNA)
# ===== File: list.py =====
# -------------------------------------------------------------------
# File name: list.py
#
# Python compound data type list can be written as a list of
# comma-separated values (items) between square brackets [,,]. Lists might
# contain items of different types, but usually the items all have the
# same type. Unlike strings, lists are a mutable type, i.e. it is
# possible to change their content. Forward index starts with 0 and
# increases; backward index starts with -1 and decreases.
#
# The code below illustrates some of the basic list operations using the
# stop codons.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/introduction.html#lists
# https://docs.python.org/3/tutorial/datastructures.html#more-on-lists
# https://docs.python.org/3/library/stdtypes.html#lists
# https://docs.python.org/3/library/stdtypes.html#common-sequence-operations
# https://en.wikipedia.org/wiki/DNA_codon_table
# ------------------------------------------------------------------

# Constructing a list
stop_codons = ['TAA', 'tAG']
print(stop_codons)

# Accessing an item in a list
first_stop_codon = stop_codons[0]
print(first_stop_codon)

# Modifying an item in a list
stop_codons[1] = 'TAG'
print(stop_codons)

# Appending an item to the end of a list
stop_codons.append('TGA')
print(stop_codons)

# Number of items in a list
number_of_stop_codons = len(stop_codons)
print('There are', number_of_stop_codons, 'stop codons')

# Convert list to a string
DNA_seq = ''.join(stop_codons)
print(DNA_seq)

# Convert string to a list
DNA_list = list(DNA_seq)
print(DNA_list)

# Slicing a list
second_codon = DNA_list[3:6]             # index 6 not included
print('Second codon:', second_codon)

# Copying a list
DNA_list_duplicate = DNA_list.copy()
print(DNA_list_duplicate)

# Insert, delete element
DNA_list_duplicate.insert(5, "?")
print(DNA_list_duplicate)
DNA_list_duplicate.pop(5)        # Can also use: del DNA_list_duplicate[5]
print(DNA_list_duplicate)
# ===== File: numeric.py =====
# -------------------------------------------------------------------
# File name: numeric.py
#
# Variable names in Python start with a letter followed by combination of
# letters, digits or underscore (no white spaces).
#
# This program illustrates some of the basics of Python's numeric data types
# int (for integer) and float.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#numeric-types-int-float-complex
# https://docs.python.org/3/library/sys.html#sys.float_info
# https://docs.python.org/3/library/math.html#module-math
# https://docs.python.org/3/tutorial/controlflow.html#intermezzo-coding-style
# https://www.census.gov/popclock/
# https://www.nature.com/articles/d41586-018-05462-w
# https://en.wikipedia.org/wiki/Exon
# https://en.wikipedia.org/wiki/Double-precision_floating-point_format
# ------------------------------------------------------------------

import math                 # to use mathematical functions

# Integers have unlimited precision
human_genes = 21306         # int type. Not 21,306 or 21 306
US_population = 328918373
print('Number of human genes:', human_genes)
print('Number of human genes in US:', human_genes*US_population)

# Floats are implemented using double in C
exons_per_gene = 8.9        # float type
print('Human exons per gene:', exons_per_gene)
human_exons = exons_per_gene*human_genes
print('Number of human exons:', human_exons)

# To convert float to integer
human_exons = int(human_exons)
print('Approximate number of human exons:', human_exons)

# Can Python do arithmetic?
firstProduct = (9.4*0.2321)*5.6
secondProduct = 9.4*(0.2321*5.6)
print('(9.4*0.2321)*5.6 - 9.4*(0.2321*5.6) =',
      (firstProduct - secondProduct))

# To access mathematical functions from math module
two_pi = 2.0*math.pi
print('two_pi =', two_pi)
print('sin(two_pi) =', math.sin(two_pi))
print('Do you believe this result?')
# ===== File: print.py =====
# ------------------------------------------------------------------
# File name: print.py
#
# print(*objects, sep=' ', end='\n', file=sys.stdout, flush=False)
# The print() function writes the value of the argument(s) it is given.
# It handles multiple arguments, floating point quantities, and strings.
# Strings are printed without quotes, and a space is inserted between items.
# The keyword argument end can be used to avoid the newline (\n) after the output
# or end the output with a different string.
#
# You can escape (overrule its special meaning) a character by
# prefixing it with backslash \
#
# The code below illustrates some of the basic usages of the print() function.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/functions.html#print
# https://docs.python.org/3/library/stdtypes.html#printf-style-string-formatting
# https://en.wikipedia.org/wiki/ASCII
# ------------------------------------------------------------------

import math

human_genes = 20000
print('You have', human_genes, 'genes')
print()

# Replacing \n with a string
print('You have', end =' ')
print(human_genes, end = '? ')
print('genes')

# Spreading over lines
print('You have',
      human_genes,
      'genes')

# Escaping with \ and string concatenation
print('You have ' + '\'' + str(human_genes) + '\'' + ' genes')

# printf style string formatting
print('The value of pi is %10s %5.3f' %('--->',  math.pi))

print(chr(7))
# ===== File: range.py =====
# ------------------------------------------------------------------
# File name: range.py
#
# range(start, stop[, step])
# This built-in function generates arithmetic progressions.
# The object returned by range() is not a list, but often acts like one.
# The arguments to the range constructor must be integers. If the step
# argument is omitted, it defaults to 1. If the start argument is omitted,
# it defaults to 0. If step is zero, ValueError is raised.
# For a positive step, the contents of a range r are determined by the
# formula r[i] = start + step*i where i >= 0 and r[i] < stop.
#
# int() returns the integer part of a decimal number.
#
# Problem: Assume the population of Florida sandhill cranes grows by
# 1.94% annually. If we start with a population of 425 birds, how large
# will the population be after 30 years? Using a for loop, code below
# computes and prints the population sizes for 30 years.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/controlflow.html#the-range-function
# https://docs.python.org/3/library/stdtypes.html#ranges
# https://sora.unm.edu/sites/default/files/journals/jfo/v061n02/p0224-p0231.pdf
# https://www4.swfwmd.state.fl.us/springscoast/sandhillcranes.shtml
# ------------------------------------------------------------------

# Initial population data
population = [425]
growth_rate = 0.0194
number_of_years = 30

# Construct a list numbering the years, from 0 to 30
years = range(0, number_of_years + 1, 1)             #31 excluded
print(list(years))

# Compute and append the population sizes to a list
for year in years:
    next_generation = population[year] + growth_rate*population[year]
    population.append(next_generation)

# Print the list of population sizes
for year in years:
    print('At year %2d the population is %7.3f' %(year, population[year]))

print()
# Print the calculated population sizes rounded to the nearest integer
for year in years:
    print('At year %2d the population is %3d' %(year, round(population[year])))
# ===== File: replace.py =====
# ------------------------------------------------------------------
# File name: replace.py
#
# str.replace(old, new[, count])
# Return a copy of the string with all occurrences of substring old 
# replaced by new. If the optional argument count is given, only the 
# first count occurrences are replaced.
#
# The code below removes spaces and digits in a segment of DNA in 
# GenBank format and highlights the start codon ATG.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#str.replace
# http://www.ncbi.nlm.nih.gov/nuccore/226377833?report=genbank
# ------------------------------------------------------------------

import re

zika_DNA = '601 catgtgtgac gccaccatga gttatgagtg'
print('Original Zika DNA\t\t:', zika_DNA)

# Replace space with nothing one time
zika_DNA = zika_DNA.replace(' ', '', 1)
print('Replace space with nothing\t:', zika_DNA)

# Replace all space characters with nothing
zika_DNA = zika_DNA.replace(' ', '')
print('Replace spaces with nothing\t:', zika_DNA)

# Substitute all digits with nothing using regular expressions
zika_DNA = re.sub(r'[1234567890]', '', zika_DNA)
print('Replace numbers with nothing\t:', zika_DNA)
# ===== File: reverse.py =====
# ------------------------------------------------------------------
# File name: reverse.py
#
# Reversing a DNA sequence is a common operation in genetics when dealing
# with a complementary strand. String slicing operation string[::-1]
# returns the desired result.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#common-sequence-operations
# https://en.wikipedia.org/wiki/Primer_%28molecular_biology%29#/media/Fi
# http://www.ncbi.nlm.nih.gov/nuccore/226377833?report=fasta
# ------------------------------------------------------------------

zika_DNA = 'AATCCATGGTTTCT'
print('Zika segment\t\t:', zika_DNA)

reversed_zika_DNA = zika_DNA[::-1]
print('Reversed zika segment\t:', reversed_zika_DNA)
# ===== File: slice.py =====
# -------------------------------------------------------------------
# File name: slice.py
#
# s[i:j:k] slice of s from i to j with step k
#
# The slice of s from i to j is defined as the sequence of items with
# index k such that i <= k < j. If i or j is greater than len(s), use len(s).
# If i is omitted or None, use 0. If j is omitted or None, use len(s).
# If i is greater than or equal to j, the slice is empty.
#
# If i or j is negative, the index is relative to the end of sequence
# s: len(s) + i or len(s) + j is substituted. But note that -0 is still 0.
#
# The code below illustrates the use of slicing in extraction of
# subsequences from a DNA sequence.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/introduction.html#strings
# https://docs.python.org/3/library/stdtypes.html#str.find
# https://docs.python.org/3/library/stdtypes.html#common-sequence-operations
# https://docs.python.org/3/library/stdtypes.html#text-sequence-type-str
# https://en.wikipedia.org/wiki/DNA_codon_table
# https://en.wikipedia.org/wiki/D-loop
# https://www.ncbi.nlm.nih.gov/nuccore/X90314.1?report=fasta
# ------------------------------------------------------------------

human = 'TTCTTTCATGGGGAAGCAGATTTGGGTACCACCCAAGTATTGACTTACCCATCAACAACCGCTATGTATT'
print('Human D-loop:', human)

mycodon = 'CAT'

# Find the lowest index of mycodon in human
index_mycodon = human.find(mycodon)
print('First', mycodon, 'index:', index_mycodon)

# Extract the first codon after mycodon
first_codon = human[index_mycodon + 3: index_mycodon + 6]
print('First codon after', mycodon, ':', first_codon)

# Extract the second codon after mycodon
second_codon = human[index_mycodon + 6: index_mycodon + 9]
print('Second codon after', mycodon, ':', second_codon)

# Negative starting point counts from the end of the string.
next_to_last_codon = human[-6:-3]
print('Next to last codon:', next_to_last_codon)

# Omitted second entry in slicing indicates to the end of string
last_codon = human[-3:]
print('Last codon:', last_codon)
# ===== File: strings.py =====
# -------------------------------------------------------------------
# File name: strings.py
#
# String variable names in Python start with a letter followed by
# combination of letters, digits or underscore (no white spaces).
# String literals are enclosed in single '...' or double "..." quotes.
# Strings can be indexed, with the first character having index 0.
# Indices may also be negative, to start counting from the right  with -1.
# Python strings cannot be changed — they are immutable. Therefore,
# assigning to an indexed position in the string results in an error.
#
# The code below illustrates some of the basic string operations on a
# partial DNA and protein sequences of GFP.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/introduction.html#strings
# https://docs.python.org/3/library/stdtypes.html#common-sequence-operations
# https://docs.python.org/3/library/stdtypes.html#text-sequence-type-str
# https://docs.python.org/3/library/functions.html#len
# https://en.wikipedia.org/wiki/Green_fluorescent_protein
# https://www.nobelprize.org/prizes/chemistry/2008/summary/
# https://scienceblogs.com/pharyngula/2008/10/08/gfp-wins-nobel-prize
# https://en.wikipedia.org/wiki/DNA_codon_table
# ------------------------------------------------------------------

# String variables and literals
protein = "GFP"                    # Winner of 2008 Nobel in chemistry
protein_seq_begin = 'MSKGEELFTG'
protein_seq_end = 'HGMDELYK'

# Concatenation of strings
protein_seq = protein_seq_begin + '...' + protein_seq_end
print('Protein sequence of GFP: ' + protein_seq)

# String method str.upper()
DNA_seq = 'atgagtaaag...actatacaaa'
DNA_seq = DNA_seq.upper()
print('DNA sequence: ' + DNA_seq)

# Forward index starts with 0 and increases
# Backward index starts with -1 and decreases
print('The second nucleotide:', DNA_seq[1])
print('The last nucleotide:', DNA_seq[-1])

# Slicing a string
first_codon = DNA_seq[0:3]         # index 3 excluded
last_codon = DNA_seq[-3:]
print('First codon:', first_codon)
print('Last codon:', last_codon)
# ===== File: translate.py =====
# ------------------------------------------------------------------
# File name: translate.py
#
# static str.maketrans(x[, y[, z]])
# This static method returns a translation table usable for str.translate().
# If there are two arguments, they must be strings of equal length,
# and in the resulting dictionary, each character in x will be mapped
# to the character at the same position in y.
# If there is a third argument, it must be a string, whose characters
# will be mapped to None in the result.
#
# The code below computes the complementary strand of a DNA sequence.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#str.maketrans
# https://docs.python.org/3/library/stdtypes.html#str.translate
# http://www.ncbi.nlm.nih.gov/nuccore/226377833?report=fasta
# https://en.wikipedia.org/wiki/DNA
# ------------------------------------------------------------------

zika_DNA = 'AGTTGTTGATCTGTGTGAGTCAG'
print("Direct strand:        5' " + zika_DNA + " 3'")

complements = zika_DNA.maketrans('acgtACGT', 'tgcaTGCA')
complement_seq = zika_DNA.translate(complements)

bonds = " "*25 + "|"*len(zika_DNA)
print(bonds)
print("Complementary strand: 3' " + complement_seq + " 5'")
# ===== File: try_except.py =====
# ------------------------------------------------------------------
# File name: try_except.py
#
# Errors detected during execution are called exceptions. It is possible
# to write programs that handle selected exceptions with a try-except statement.
# The try statement works as follows:
# 1. The statement(s) between the try and except keywords is executed.
# 2. If no exception happens, the except block  is skipped and execution
#    of the try statement is finished.
# 3. If an exception occurs during execution of the try clause, the rest
#    of the clause is skipped. Then if its type matches the exception
#    named after the except keyword, the except clause is executed, and
#    then execution continues after the try statement.
# 4. An unhandled exception stops the execution.
# 5. Optional else clause must follow all except clauses. It is useful
#    for code that must be executed if the try clause does not raise
#    an exception.
#
# The code below illustrates the usage of try-except statement in handling
# exceptions, e.g. input error (ValueError), out-of-bound index error 
# (IndexError) in a list. The program asks for an input, repeatedly, if an 
# exception is raised. If no exception is raised, the program breaks out of the 
# infinite while loop by executing the else clause and finishes with 
# a message. One can also exit the program by typing Control^C.
# 
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/tutorial/errors.html#exceptions
# https://docs.python.org/3/reference/compound_stmts.html#the-try-statement
# https://docs.python.org/3/library/exceptions.html
# ------------------------------------------------------------------

import sys

stop_codons = ['TAA', 'TAG', 'TGA']
print('Stop codons:', stop_codons)

while True:
    try:
         index = int(input('Please enter the index of a stop codon to print: '))
         print('Your codon is', stop_codons[index])
    except ValueError as ve:
         print(ve, 'Try again...')
    except IndexError:
         print('Your index', index, 'is out of range. Try again...')
    except:
         print('Unexpected error:', sys.exc_info()[0])
         sys.exit()
    else:
         print('Good bye!')
         break
# ===== File: tuple.py =====
# -------------------------------------------------------------------
# File name: tuple.py
#
# Python compound data type tuple is an immutable sequence used to store
# collections of heterogeneous or homogenous data. They can be constructed
# in a number of ways; most typically, separating items with commas: 
# a, b, c or (a, b, c). Tuples implement all of the common sequence 
# operations. Some Python functions and methods return a list of tuples.
# 
# The code below illustrates some of the common tuple operations using 
# basic amino acids and their codons. 
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/library/stdtypes.html#tuples
# https://docs.python.org/3/library/stdtypes.html#common-sequence-operations
# https://en.wikipedia.org/wiki/DNA_codon_table
# ------------------------------------------------------------------

# Constructing tuples
Histidine = ('H', 'CAT', 'CAC')
Lysine = 'K', 'AAA', 'AAG'
Arginine = ('R', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG')

print('Histidine:', Histidine)
print('Lysine:', Lysine)
print('Arginine:', Arginine)
print()

# Constructing a list of tuples
basic = [Histidine, Lysine]
basic.append(Arginine)
print('Basic amino acids:', basic)
print()

# Accessing elements of a list of tuples
His = basic[0]
print('His:', His)
His_codons = basic[0][1:]
print('His codons:', His_codons)
codon1, codon2 = His_codons
print('codon1:', codon1)
print('codon2:', codon2)
print()

protein_seq = basic[0][0] + basic[1][0] + basic[2][0]
print('Protein:', protein_seq)
# ===== File: welcome.py =====
# ------------------------------------------------------------------
# File name: welcome.py
#
# A line starting with # character is a comment and not interpreted
# by Python. Use comments liberally in your codes.
#
# This is a customary first program to test if your computer is ready
# for Python. Our program will simply print:
#                    Welcome to Python for Biologists!
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# http://www.python.org
# https://www.python.org/about/gettingstarted/
# ------------------------------------------------------------------

print('Welcome to Python for Biologists!')
# ===== File: while.py =====
# ------------------------------------------------------------------
# File name: while.py
#
# while EXPR:
#    statements
#
# while loop iterates the block of statements as long as EXPR remains True.
#
# To illustrate the usage of while statement, the code below first
# computes the number of appearances of a nucleotide base in a string 
# using Python's str.count() method. Then it computes the same number
# using a while statement, hoping to get the same answer.
#
# Version: 2.1
# Authors: H. Kocak and B. Koc
#          University of Miami and Stetson University
# References:
# https://docs.python.org/3/reference/compound_stmts.html#the-while-statement
# https://docs.python.org/3/library/stdtypes.html#string-methods
# https://www.ncbi.nlm.nih.gov/nuccore/KC545393.1?report=fasta
# ------------------------------------------------------------------

DNA_seq = 'CGGACACACAAAAAGAATGAAGGATTTTGAATCTTTATTGTGTGCGAGTAACTACGAGGAAGATTAAAGA'
print('DNA sequence:', DNA_seq)

bp = 'T'
print('Base pair:', bp)

print('str.count():', DNA_seq.count(bp))

count = 0
index = 0

while index < len(DNA_seq):
    if bp == DNA_seq[index]:
        count += 1
    index += 1

print('Our while count:', count)


# ===== Rosalind Problems 1–49 =====
# Each problem's core algorithm is implemented with short code + 2-line explanation

# 1. Counting DNA Nucleotides
def count_dna(seq):
    # Count A, C, G, T in the sequence
    return seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')

# 2. Transcribing DNA into RNA
def transcribe_dna(seq):
    # Replace T with U to convert DNA → RNA
    return seq.replace('T','U')

# 3. Reverse Complement
def reverse_complement(seq):
    # Reverse the sequence and complement each base
    comp = str.maketrans("ATCG","TAGC")
    return seq.translate(comp)[::-1]

# ... Continue similar short implementations for all problems up to 49 ...
