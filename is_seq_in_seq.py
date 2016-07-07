#!/bin/bash

# imports
import os, tkFileDialog, datetime, subprocess, csv
from Tkinter import Tk


"""
This script aligns input from a designated input CSV file in the form of labeled nucleotide
sequences, and compares them to a folder (or directory tree) or sequencing results
to try to identify matches between the designated input sequences and the queried
sequences.

The format for the input CSV is as follows:

+----------------+----------------+
| Sequence Label | Input Sequence |
+----------------+----------------+
|     seq 1      | CTACGTGAGCGATG |
|     seq 2      | TCTGGCGCTCGACA |
|     seq 3      | GGCGCTAATATGCG |
|      ...       |       ...      |
|     seq N      | NNNNNNNNNNNNNN |
+----------------+----------------+

The script uses simple string regex'ing, but could be modified to use sequence alignment
algorithms such as the Clustal Omega cli. 

"""

class SeqSearch:
	def __init__(self):
		self.initial_timestamp = self.get_timestamp()

		# File Paths
		self.root_directory = "/Users/sdodd/Documents/Data/SeqSearch/"
		self.run_directory = self.root_directory+"{}_seq_search/".format(self.initial_timestamp)
		self.input_csv_path = "/Users/sdodd/Desktop/psm.csv"#str(raw_input("Input CSV location: "))
		self.target_sequence_path = "/Users/sdodd/Documents/Data/Sequencing/Results/2016-07-01/"#str(raw_input("Target Sequences Folder: "))

		# Search data
		"""
		Hack work-around because I don't have time.
		Fix this Fix this Fix this Fix this Fix this Fix this Fix this Fix this 
		
		self.reverse_primer = True # TODO TODO REMOVE THIS AND READ THIS VALUE PER SEQUENCE FROM CSV
		
		Fix this Fix this Fix this Fix this Fix this Fix this Fix this Fix this 
		"""

		self.identity_cutoff = 40
		self.searching = True
		self.minimum_search_length = 10
		self.input_sequences = {}
		self.target_sequences = {}
		self.match_found = False
		self.matched_sequences = []

		# makes all of the directories
		self.makedir(self.root_directory)
		self.makedir(self.run_directory)

		# output data
		self.log = ""
		# self.successes = "" 										#unused
		# self.failures = "" 										#unused


	"""
	Returns a timestamp in the form of year_month_day_hour_minute_second
	"""
	def get_timestamp(self):
		today = datetime.datetime.today().strftime('%y_%m_%d')
		hour = datetime.datetime.today().time().hour
		minute = datetime.datetime.today().time().minute
		second = datetime.datetime.today().time().second

		timestamp = "{}_{}_{}_{}".format(today, hour, minute, second)
		return timestamp

	"""
	Makes a directory if it does not already exist
	"""
	def makedir(self, dir_to_make):
		if not os.path.exists(dir_to_make):
			os.makedirs(dir_to_make)
		else:
			print "directory already exists: {}".format(dir_to_make)

	"""
	Adds string to the log
	"""
	def log_line(self, string_to_log):
		self.log += str(string_to_log) + "\n"

	"""
	------------------------------------------------------------------------
	------------------------------------------------------------------------
	---------------------------- MEAT METHODS ------------------------------
	------------------------------------------------------------------------
	------------------------------------------------------------------------
	"""

	"""
	This method uses location of the input csv file, along with the location of the 
	sequence folder with the sequences to align to, and creates an unaligned FASTA 
	file to be used by the Clustal Omega algorithm for alignment.
	"""
	def get_input_sequences(self):
		# read the data from the input file and save the sequence name / seq
		with open(self.input_csv_path, 'rb') as input_file:
			reader = csv.reader(input_file, delimiter=',')
			for row in reader:
				if not row[0] == "Label":
					self.input_sequences[row[0]] = row[1]
					print "adding sequence {} | {} to input".format(row[0], row[1])
			input_file.close()

	"""
	This method searches in the target sequence for the input sequence. If the
	whole input sequence is not found, the input sequence is chopped into
	contiguous subsections (>=10bp), which are then searched for in the target
	sequence. 
	"""
	def search_for_sequence(self):

		for target_sequence in self.target_sequences:

			print "target: {}-------------------------------------".format(target_sequence)

			target_fasta = self.target_sequences[target_sequence]
		
			if len(target_fasta) == 0:
				raise ValueError("Target Sequence is empty")
			else:
				for seq in self.input_sequences:
					print "input sequence: {} | {}".format(seq, self.input_sequences[seq])
					self.match_found = False
					fasta_sequence = self.input_sequences[seq]
					if len(fasta_sequence) < self.minimum_search_length:
						raise ValueError("seq length must be > 10bp")
					
					# search while we haven't found a match
					while not self.match_found:

						# get all of the sequence contigs in all positions of size
						# from the full sequence to 10bp
						for x in range(0, len(fasta_sequence) - self.minimum_search_length):
							length_of_current_search_string = len(fasta_sequence) - x
							contigs_to_search_for = []
							for y in range(0, len(fasta_sequence) - length_of_current_search_string):
								current_contig = fasta_sequence[y:length_of_current_search_string + y]
								if len(current_contig) >= self.minimum_search_length:
									contigs_to_search_for.append(current_contig)

							# search for the search seqs in the target
							for search_contig in contigs_to_search_for:
								if search_contig in target_fasta:
									print "match found: {} | {}".format(seq, search_contig)
									match_result = SequenceMatchResult(seq, search_contig, target_sequence)
									self.matched_sequences.append(match_result)
									self.match_found = True
									break
								else:
									pass

						# if no match is found for this input sequence in the given target,
						# we set the searching variable to false so we can move on to the
						# next input sequence
						self.match_found = True

	"""
	Adds all of the target sequences in the selected file folder path to the
	search.
	"""
	def get_target_sequences(self, target_sequence_path):

		# if the path is a directory of sequences
		if os.path.isdir(target_sequence_path):

			print "list dir: {}".format(os.listdir(target_sequence_path))
			for file_path in os.listdir(target_sequence_path):

				print "file path: {}".format(file_path)

				# if the path is another nested directory
				if os.path.isdir(target_sequence_path+file_path):
					print "getting sequences from nested dir: {}".format(target_sequence_path+file_path)
					self.get_target_sequences(target_sequence_path+file_path)

				if file_path[-4:] == ".seq":
					# format the file path
					if target_sequence_path[-1] == "/":
						full_file_path = target_sequence_path + file_path
						print 'adding target file: {}'.format(full_file_path)
						self.get_target_sequence(full_file_path)
					else:
						full_file_path = target_sequence_path + "/" + file_path
						print 'adding target file: {}'.format(full_file_path)
						self.get_target_sequence(full_file_path)
		# if the path is a single sequence file
		else:
			self.get_target_sequence(self.target_sequence_path)


	"""
	Adds a value to self.target_sequences {seq name : sequence} for the target sequence. The
	data for this dict is parsed from the target file that is input at the creation of
	the alignment object.
	"""
	def get_target_sequence(self, target_file_path):

		target_name = "target"
		target_sequence = ""

		with open(target_file_path, "rb") as target_file:
			target_data = target_file.read()
			if ">" in target_data:
				target_name = target_data.splitlines()[0]
				target_sequence = ""
				for line in target_data.splitlines()[1:]:
					target_sequence+=line

				self.target_sequences[target_name] = target_sequence
				print self.target_sequences


			else:
				target_name = os.path.split(target_file_path)[-1]
				target_sequence = target_data.replace("\n","")

				self.target_sequences[target_name] = target_sequence
				print self.target_sequences

	"""
	Formats the results for easy output
	"""


	"""
	Returns the results for the search in the form of 
	"""
	def return_results(self):

		selected_results = []

		self.log_line("PERCENT IDENTITY CUTOFF: {}".format(self.identity_cutoff))

		for index, result in enumerate(self.matched_sequences):
			search_fasta = self.input_sequences[result.search_seq]
			target_fasta = self.target_sequences[result.target_seq]
			percent_hit_identity = float(len(result.match_contig)) / len(search_fasta) * 100

			if percent_hit_identity > self.identity_cutoff and result.search_seq not in selected_results:

				selected_results.append(result.search_seq)

				self.log_line("--------RESULT {} ----------------------------".format(len(selected_results)))
				self.log_line("Search sequence: {}".format(result.search_seq))
				self.log_line("Search FASTA: {}".format(search_fasta))
				self.log_line("Match FASTA: {}".format(result.match_contig))
				self.log_line("Target sequence: {}".format(result.target_seq))
				self.log_line("Target FASTA: {}".format(target_fasta))
				self.log_line("Percent Hit Identity: {}".format(percent_hit_identity))
				self.log_line("----------------------------------------------")

		print self.log

		# output log to results file
		with open(self.run_directory+"results.txt", "wb") as save_file:
			save_file.write(self.log)
			save_file.close()



	"""
	Runs the show
	"""
	def run(self):
		self.get_target_sequences(self.target_sequence_path)
		self.get_input_sequences()
		self.search_for_sequence()
		self.return_results()



"""
Class to hold the result of a contig match from a sequence search
"""

class SequenceMatchResult:
	def __init__(self, search_seq, match_contig, target_seq):
		self.search_seq = search_seq
		self.match_contig = match_contig
		self.target_seq = target_seq





def main():
	alignment = SeqSearch()
	alignment.run()

if __name__ == "__main__":
	main()
















