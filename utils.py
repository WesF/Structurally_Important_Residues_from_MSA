import os
import os.path
import logging
import math
import xml.dom.minidom
from Bio.SeqIO.PdbIO import PdbSeqresIterator

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Align.Applications import ClustalwCommandline
from Bio.PDB import PDBParser, PPBuilder

class Step(object):
	name = ''
	output_files = [] # override with the output files
	
	def clean(self, state):
		for f in self.output_files:
			path_to_delete = state[f]
			if os.path.exists(path_to_delete):
				logging.debug('deleting %s' % (path_to_delete))
				os.remove(path_to_delete)

	def run(self, state, input):
		raise NotImplementedError()

def dump_matrix_to_file(matrix, path):
	# save to file
	with open(path, "w") as save_file:
		for row in matrix:
			line = "\t".join(["%8f" % (float(x)) for x in row])
			save_file.write(line)
			save_file.write("\n")

		
class Step1(Step):
	name = 'Step1'
	output_files = ['distances_matrix', 'reference_sequence']
	
	def run(self, state):
		parser = PDBParser()
		with file(state.input_pdb,"r") as pdbfile:
			x = PdbSeqresIterator(pdbfile)
			pdbsequence = list(x)[0]
			with file(state.reference_sequence, "w") as ref_seq_file:
				ref_seq_file.write('>reference\n')
				for c in pdbsequence:
					ref_seq_file.write(str(c))
				logging.info('sequence was extracted from pdb file. sequence length=%d' % (len(pdbsequence)))
		
		structure = parser.get_structure('myPDB',state.input_pdb)
		
		
		first_model= structure.__iter__().next()
		first_chain = first_model.__iter__().next()

		max_position = 0
		for residue in first_chain:
			try:
				max_position = max(max_position, residue.id[1])
			except KeyError:
				pass
		
		positions = [0] * max_position
		print "length is"+str(max_position)
		
		i=0
		for residue in first_chain:
			try:
				r1 = residue['CA']
				print residue.id[1]
				#positions[residue.id[1]-1] = r1
				positions[i] = r1
				i+=1
			except KeyError:
				pass
		distance_matrix=[]
		for residue_1 in positions:
			foo = []
			for residue_2 in positions:
				if residue_1 == 0 or residue_2 == 0:
					foo.append(0)
				else:
					foo.append(residue_1 - residue_2)
			distance_matrix.append(foo)
		dump_matrix_to_file(distance_matrix, state.distances_matrix)

class Step2(Step):
	name = 'manual blast'
	output_files = ['blast3fasta']
	
	def run(self, state):
		print 'Manual blast'
		print '============'
		print 'Input : %s' % (state.reference_sequence)
		print 'Output: %s' % (state.blast3fasta)
		print 'Press enter to continue'
		
		# add the reference to the result		 
		raw_input()
		with file(state.blast3fasta, "a") as output:
		  with file(state.reference_sequence, "r") as ref_seq_file:
			output.write(ref_seq_file.read())

class MsaStep(Step):
	input = None
	output = None
	format = 'CLUSTAL'
	def run(self, state):
		POSSIBLE_PATHS= ["/Users/Wes/Desktop/clustalo","c:\\Program Files (x86)\\ClustalW2\\clustalw2.exe"]
		selected_path = None
		for path in POSSIBLE_PATHS:
			if os.path.exists(path):
				selected_path = path
				break		 
		cline = ClustalwCommandline(selected_path, infile=state[self.input], outfile=state[self.output], output=self.format)
		cline()
		

class Step3(MsaStep):
	name = 'Step3'
	output_files = ['msa']
	input = 'blast3fasta'
	output = output_files[0]
	format = 'FASTA'

class Step4(Step):
	name = 'Step4'
	output_files = ['correlations']
		
	def run(self, state):
		with file(state.msa) as input_file:
			sequences = list(SeqIO.parse(input_file, format="fasta"))
		ref = [s for s in sequences if s.name=='reference']
		non_ref = [s for s in sequences if s.name != 'reference']
		sequences = ref
		sequences.extend(non_ref)
		
		# Translation
		# A - Positive
		# B - Negative
		# C - Polar
		# D - Nonpolar
		# E - Aromatic
		# F - Proline
		# X - Unidentified
		# - - Gap
		DICT = {
		  'A' : 'D',   # Alanine
		  'R' : 'A',   # Arginine
		  'B' : 'C',   # Aspartic Acid or Aspargine
		  'N' : 'C',   # Aspargine
		  'D' : 'B',   # Aspartic Acid
		  'C' : 'C',   # Cysteine
		  'Q' : 'C',   # Glutamine
		  'E' : 'B',   # Glutamic Acid
		  'Z' : 'C',   # Glutamic Acid or Glutamine
		  'G' : 'D',   # Glycine
		  'H' : 'A',   # Histidine
		  'I' : 'D',   # Isoleucine
		  'L' : 'D' ,  # Leucine
		  'K' : 'A',   # Lysine
		  'M' : 'D',   # Methionine
		  'F' : 'E',   # Phenylalanine
		  'P' : 'F',   # Proline
		  'S' : 'C',   # Serine
		  'T' : 'C',   # Threonine
		  'W' : 'E',   # Tryptophan
		  'Y' : 'E',   # Thyrosine
		  'V' : 'D',   # Valine
		  'X' : 'X',   # Unidentified
		  '-' : '-',   # Gap
		}
		
		for r in sequences:
			for i in range(len(r)):
				try:
				  y = DICT[r[i]]
				except:
				  pass
				  print r
				  print i
				  print r[i]
				  print ''
		#sequences = [[DICT[i] for i in row] for row in sequences]
		
		len_sequences = len(sequences)
		N = len(sequences[0])
		logging.info('sequence length is %d' % (N))
		logging.info('found %d sequences' % (len_sequences))
		self.calc_factorial_dict(len_sequences)
		
		# create 2-dim array N x N
		correlations = [[0 for x in xrange(N)] for x in xrange(N)]
		
		calcs = N * (N+1) / 2
		current_calc = 0
		last_precent = 0
		# go over positions pairs
		for pos_i in xrange(N):
			for pos_j in xrange(pos_i): # j run from 0 to pos_i-1
				correlations[pos_i][pos_j] = self.calc_correlation(sequences, len_sequences, pos_i, pos_j)
				current_calc+=1
				new_precent = current_calc * 100 / calcs
				if new_precent != last_precent:
					logging.info('%d%% calculated' % (new_precent))
					last_precent = new_precent
		
		
		# save to file
		dump_matrix_to_file(correlations,state.correlations)  
	
	def calc_factorial_dict(self, N):
		self.factorial = {0:1}
		value = 1
		for i in xrange(1,N+1):
			value *= i
			self.factorial[i] = value
	
	def calc_correlation(self, sequences, len_sequences, pos_i, pos_j):
		# assuming the first sequence is the reference sequence
		ref_i = sequences[0][pos_i]
		ref_j = sequences[0][pos_j]
		
		na = 0
		nb = 0
		nba = 0
		for nref in xrange(1,len_sequences):
			nref_i = sequences[nref][pos_i]
			nref_j = sequences[nref][pos_j]
			if ref_i == nref_i and ref_j == nref_j:
				na+=1
				nb+=1
				nba+=1
			elif ref_i == nref_i and ref_j != nref_j:
				na+=1
			elif ref_i != nref_i and ref_j == nref_j:
				nb+=1
		
		expected_nba = float(na) * nb / len_sequences
		if expected_nba < nba:
			min_index = nba
			max_index = na
			multiply_by = -1
		elif expected_nba > nba:
			min_index = 0
			max_index = nba
			multiply_by = 1
		else:
			return 0
		
		sum = long(0)
		for i in xrange(min_index, max_index+1):
			first_part = long(self.factorial[na] / (self.factorial[i] * self.factorial[na-i]))
			second_part = (long(nb)/len_sequences) ** i
			third_part = (1-long(nb)/len_sequences) ** (na-i)
			
			sum += first_part * second_part * third_part
		return multiply_by * math.log10(sum)

class Step5(Step):
	name = 'Step5'
	output_files = ['correlations_reduced']

	def run(self, state):
		with file(state.msa) as input_file:
			sequences = list(SeqIO.parse(input_file, format="fasta"))
		ref = [s for s in sequences if s.name=='reference']
		
		if len(ref) != 1:
			logging.error('number of refernce sequences != 1')
			#print ref[0] 
			return
		ref = ref[0]

		correlations = []
		with file(state.correlations) as input_file:
			for line in input_file:
				correlations.append(line.split("\t"))
		ref = ''.join(list(ref))
		while ref.find('-') != -1:
			index_to_remove = ref.index('-')
			correlations = correlations[:index_to_remove] + correlations[index_to_remove+1:]
			correlations = [row[:index_to_remove] + row[index_to_remove+1:] for row in correlations]
			ref = ref[:index_to_remove] + ref[index_to_remove+1:]
		dump_matrix_to_file(correlations,state.correlations_reduced)
		
	
class Step6(Step):
	name = 'Step6'
	output_files = ['correlation_for_pair']
	
	def run(self, state):
		correlations = []
		cells = []
		with file(state.correlations_reduced) as input_file:
			for line in input_file:
				correlations.append(line.split("\t"))
		for pos_i in xrange(len(correlations)):
			row_i = correlations[pos_i]
			for pos_j in xrange(len(row_i)):
				val = row_i[pos_j]
				cells.append((float(val), pos_i, pos_j))
		cells.sort()
		# save to file
		with open(state.correlation_for_pair, "w") as save_file:
			for row in cells:
				save_file.write(str(row))
				save_file.write("\n")
				
class Step7(Step):
	name = 'Step7'
	output_files = []	 
	def run(self, state):
		correlations = []
		with file(state.correlations_reduced) as input_file:
			for line in input_file:
				correlations.append(line.split("\t"))
		distances = []
		with file(state.distances_matrix) as input_file:
			for line in input_file:
				distances.append(line.split("\t"))		  
		if (len(correlations) != len(distances)):
			logging.error('correlations has %d lines. distances has %d lines' % (len(correlations), len(distances)))
			return #cant continue
		else:
			logging.info('correlations has %d lines. distances has %d lines' % (len(correlations), len(distances)))
			
		x = []
		y = []
		c = []
		
		# Translation
		# A - Positive
		# B - Negative
		# C - Polar
		# D - Nonpolar
		# E - Aromatic
		# F - Proline
		# X - Unidentified
		# - - Gap
		DICT = {
		  'A' : 'D',   # Alanine
		  'R' : 'A',   # Arginine
		  'B' : 'C',   # Aspartic Acid or Aspargine
		  'N' : 'C',   # Aspargine
		  'D' : 'B',   # Aspartic Acid
		  'C' : 'C',   # Cysteine
		  'Q' : 'C',   # Glutamine
		  'E' : 'B',   # Glutamic Acid
		  'Z' : 'C',   # Glutamic Acid or Glutamine
		  'G' : 'D',   # Glycine
		  'H' : 'A',   # Histidine
		  'I' : 'D',   # Isoleucine
		  'L' : 'D' ,  # Leucine
		  'K' : 'A',   # Lysine
		  'M' : 'D',   # Methionine
		  'F' : 'E',   # Phenylalanine
		  'P' : 'F',   # Proline
		  'S' : 'C',   # Serine
		  'T' : 'C',   # Threonine
		  'W' : 'E',   # Tryptophan
		  'Y' : 'E',   # Thyrosine
		  'V' : 'D',   # Valine
		  'X' : 'X',   # Unidentified
		  '-' : '-',   # Gap
		}
		
		PAIRS_TO_COLOR = {
			('A','A') : 'r',
			('B','B') : 'r',
			('A','B') : 'b',
			('B','A') : 'b',
			('D','D') : 'g',
			
		}
		
		with file(state.reference_sequence) as refseq:
			sequence = list(SeqIO.read(refseq, format="fasta"))
			
			
		for i in xrange(len(correlations)):
			#check to make sure length of sequence is same as length of correlation matrix
			resA=sequence[i]
			for j in xrange(i):
				#select high correlation pairs only
				#if float(correlations[i][j])<40:
					#continue
				if i == 28 or i== 63 or i== 67 or i== 42 or i== 106:
					x.append(float (correlations[i][j]))
					y.append(float(distances[i][j]))
					resB=sequence[j]
					resA_trans = DICT[resA]
					resB_trans = DICT[resB]
					color_to_append= PAIRS_TO_COLOR.get((resA_trans,resB_trans),'k')
					if i == 28:
						c.append('m')
					elif i == 42:
						c.append('y')
						
					elif i == 63:
						c.append('r')
					elif i == 67:
						c.append('b')
					elif i == 106:
						c.append('g')
					else:
						c.append('w')
				
				#c.append(color_to_append)
				
		print len(x)
		print"\n"
		print len(y)
		#try:
		plt.scatter(x, y,c=c)
		#except ValueError:
			#print 'xxx' 
	
		plt.show()

STEPS = [Step1, Step2, Step3, Step4, Step5, Step6, Step7]

##### IGNORE BEYOND THIS POINT
# class Step2(Step):
	# name = 'psiblast(1)'
	# output_files = ['blast1xml']
	
	# def run(self, state):
		# input = state.reference_sequence
		# cline = NcbipsiblastCommandline(remote=True, query=state.reference_sequence, db='nr', outfmt=5, out=state.blast1xml, max_target_seqs=50)
		# cline()

# class Xml2FastaStep(Step):
	# input = None
	# Output = None

	# def run(self, state):
		# doc = xml.dom.minidom.parse(state[self.input])
		# with file(state[self.output], "w") as save_file:
			# i = 1
			# for element in doc.getElementsByTagName('Hsp_hseq'):
				# save_file.write('>%d\n%s\n' % (i,str(element.childNodes[0].nodeValue)))
				# i += 1

# class Step3(Xml2FastaStep):
	# name = 'xml2fasta(1)'
	# output_files = ['blast1fasta']
	# input = 'blast1xml'
	# output = output_files[0]
		
# class Step2Old(Step):
	# name = 'Step2'
	# output_files = ['blast']
	
	# def run(self, state):
		# input = state.reference_sequence
		# logging.info('running BLAST for %s' % (input))
		# record = SeqIO.read(file(input), format="fasta")
		# result_handle = NCBIWWW.qblast("blastp", "nr", record.seq, hitlist_size=500)
		
		# save_file = open(state.blast, "w")
		# # put the reference
		# save_file.write('>reference\n%s\n' % (str(record.seq)))
		
		# blast_records = NCBIXML.parse(result_handle)
		
		# alignment_number = 1
		# for i,blast_record in enumerate(blast_records):
			# for alignment in blast_record.alignments:
				# for hsp in alignment.hsps:
					# save_file.write('>%d\n%s\n' % (alignment_number, hsp.query))
					# alignment_number += 1
		# save_file.close()
		# result_handle.close()


# class Step4(MsaStep):
	# name = 'fasta2aln(1)'
	# output_files = ['blast1aln']
	# input = 'blast1fasta'
	# output = output_files[0]

# class PsiBlastIterationStep(Step):
	# input = None
	# output = None
	
	# def run(self, state):
		# cline = NcbipsiblastCommandline(remote=True, in_msa=state[self.input], db='nr', outfmt=5, out=state[self.output], max_target_seqs=50)
		# cline()

# class Step5(PsiBlastIterationStep):
	# name = 'psiblast(2)'
	# input = 'blast1aln'
	# output_files = ['blast2xml']
	# output = output_files[0]
	
# class Step6(Xml2FastaStep):
	# name = 'xml2fasta(2)'
	# output_files = ['blast2fasta']
	# input = 'blast2xml'
	# output = output_files[0]
	
# class Step7(MsaStep):
	# name = 'fasta2aln(2)'
	# output_files = ['blast2aln']
	# input = 'blast2fasta'
	# output = output_files[0]

# class Step8(PsiBlastIterationStep):
	# name = 'psiblast(3)'
	# input = 'blast2aln'
	# output_files = ['blast3xml']
	# output = output_files[0]
	
# class Step9(Xml2FastaStep):
	# name = 'xml2fasta(3)'
	# output_files = ['blast3fasta']
	# input = 'blast3xml'
	# output = output_files[0]
