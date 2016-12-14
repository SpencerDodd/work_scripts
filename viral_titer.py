"""
This program calculates the concentration of a viral titer based on the results of the Broad 
Institute's Genetic Perturbation Platform's Lentiviral Titration Protocol.
"""
import xlrd
import numpy
import numpy.polynomial.polynomial
import matplotlib.pyplot as plt
import sys
import os
import Tkinter, tkFileDialog
import datetime
from scipy.optimize import curve_fit

"""
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Variables ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
"""

# Constants

# Initialization


"""
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Objects /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
"""
class Plate:

	"""
	Result workbook coordinates for the top left value of read values
		- Plate is 12 values wide and 8 values high
			i.e. x size = 12, y size = 8
	"""
	x = 2
	y = 41

	def __init__(self, file_path):
		self.file_path = file_path
		self.column_values = self.get_column_values()
		self.starting_concentration = self.get_starting_concentration()
		self.standard_titration = self.create_standard_titration()
		self.standard_curve = self.create_standard_curve()
		self.plot = self.plot_standards()
		self.plot_exponential_fit = None

	def __init__(self, file_path, starting_concentration):
		self.file_path = file_path
		self.column_values = self.get_column_values()
		self.starting_concentration = starting_concentration
		self.standard_titration = self.create_standard_titration()
		self.standard_curve = self.create_standard_curve()
		self.plot = self.plot_standards()
		self.plot_exponential_fit = None

	"""

	Reads values from an Excel .xlsx file. 
	The method places the collected values in a nested array structure.

		+----+----+----+----+----+----+----+----+----+----+----+----+
		| C  | C  | C  | C  | C  | C  | C  | C  | C  | C  | C  | C  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| o  | o  | o  | o  | o  | o  | o  | o  | o  | o  | o  | o  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| l  | l  | l  | l  | l  | l  | l  | l  | l  | l  | l  | l  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| u  | u  | u  | u  | u  | u  | u  | u  | u  | u  | u  | u  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| m  | m  | m  | m  | m  | m  | m  | m  | m  | m  | m  | m  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| n  | n  | n  | n  | n  | n  | n  | n  | n  | n  | n  | n  |
		+----+----+----+----+----+----+----+----+----+----+----+----+
		| 01 | 02 | 03 | 04 | 05 | 06 | 07 | 08 | 09 | 10 | 11 | 12 |
		+----+----+----+----+----+----+----+----+----+----+----+----+

		where column_values[0] == Column01 (above)

	The values are stored in the Plate field 'self.column_values'

	"""
	def get_column_values(self):
		wb = xlrd.open_workbook(self.file_path)
		s = wb.sheet_by_index(0)

		columns = []

		for i in range(Plate.x, Plate.x+12): 			# length + 1 due to non-inclusive upper bounds
			current_column = []
			
			for j in range(Plate.y, Plate.y+8):			# length + 1 due to non-inclusive upper bounds
				cell_value = s.cell(j, i).value

				"""
				Parse out 'OVERFLOW' string and other non-float or broken data here
				"""
				filtered_cell_value = cell_value

				current_column.append(filtered_cell_value)

			columns.append(current_column)

		return columns

	"""
	Collects input for the starting concentration of the standard curve virus
	"""
	def get_starting_concentration(self):
		starting_concentration = None
		while starting_concentration is None:
			try:
				user_input = float(input("What was the starting pRosetta Concentration? (pfu/mL): "))

			except KeyboardInterrupt:
				print ("\nExiting ...")
				sys.exit()
			
			except:
				print ("Please enter integers only to represent the concentration")

			else:
				starting_concentration = user_input

		return starting_concentration

	"""
	Creates an array of the standard curve concentration titration values given the starting
	concentration of the control virus pRosetta
	"""
	def create_standard_titration(self):
		standard_titration = []
		for iteration in range(0, 8):
			current_concentration = self.starting_concentration
			
			for x in range(0, iteration):
				current_concentration = current_concentration / 2

			standard_titration.append(current_concentration)

		print standard_titration
		return standard_titration

	"""
	Creates the standard curve by finding the average value for each column that contains
	standard virus. Returns a dict with the average standard curve in the following format:
		dict[dilution] = (fluorescence value)
	"""
	def create_standard_curve(self):
		number_of_standard_curve_columns = int(input("How many columns are pRosetta?: "))
		standard_columns = []
		for i in range(0, number_of_standard_curve_columns):
			column_to_add = self.column_values[int(input("Index of pRosetta column (0-based): "))]
			standard_columns.append(column_to_add)

		standard_averages = []

		# calculate the average value for each replicate of a standard's concentration
		for i in range(0, 8):
			sum_for_concentration = 0

			for column in standard_columns:
				sum_for_concentration += column[i]

			standard_averages.append(sum_for_concentration / len(standard_columns))

		standard_curve = {}

		# create a dict of the standard average values
		for index, average in enumerate(standard_averages):
			standard_curve[self.standard_titration[index]] = average

		return standard_curve

	"""
	Creates a plot of the standard curve
	Sets the value of self.plot_exponential_fit
	Sets the value of self.back_calc_exponential_fit
	"""
	def plot_standards(self):

		# create arrays that contain tuples of the x and y values for the standard curve
		x_vals = self.standard_titration
		y_vals = []

		# y vals are in order of x_vals (highest concentration to lowest)
		for concentration in x_vals:
			y_vals.append(self.standard_curve[concentration])

		# reverse the arrays so that they are in the correct order (i.e. low -> high)
		x_vals.reverse()
		y_vals.reverse()

		# create numpy arrays from our value arrays
		np_x_vals = numpy.array(x_vals)
		np_y_vals = numpy.array(y_vals)
	 
		# gives us the equation 
		self.plot_exponential_fit = get_exponential_fit(np_x_vals, np_y_vals, f_use="plotting")
		self.back_calc_exponential_fit = get_exponential_fit(np_x_vals, np_y_vals, f_use="back-calc")

		# calculate x and y values of the curve that best fits the dataset
		x_curve = numpy.linspace(np_x_vals[0], np_x_vals[-1], 50)
		y_curve = (lambda x: numpy.exp(self.plot_exponential_fit(numpy.log(x))))(x_curve)

		# set the title of the figure
		fig = plt.figure()
		# set the title of the figure
		fig.canvas.set_window_title(self.file_path)

		# plot the curve (blue circles)
		plt.plot(np_x_vals, np_y_vals, 'o', x_curve, y_curve)

		# label axes
		plt.ylabel("Fluorescence (RFU)")
		plt.xlabel("Viral Concentration (PFU/mL)")

		# return the created plot
		return plt

	"""
	Plots the values of the unknowns
	"""
	def plot_unknowns_on_standard_curve(self):
		
		for column in self.column_values:
			ordered_column = column[::-1]
			y_curve = numpy.array(ordered_column)
			x_curve = (lambda x: numpy.exp(self.back_calc_exponential_fit(numpy.log(x))))(y_curve)

			if self.plot is not None:
				print (x_curve)
				print (y_curve)

				self.plot.plot(x_curve, y_curve)
			
			else:
				print ("Plot needs to exist with standards before we add unknowns.")

	"""
	Plots the unknowns with their own curves on the same plot as the standard curve. Assumes
	the x values of the unknowns (concentration) to be the same as the standard
	"""
	def plot_unknowns_with_standard_curve(self):
		
		for column in self.column_values:
			# get the values for the unknown reads in ascending order
			ordered_unknown_reads = column[::-1]
			# get the values for the standard curve in ascending order
			ordered_standard_curve = self.standard_titration
			ordered_standard_curve.reverse()

			# turn the values for unknown into numpy arrays for plotting
			x_curve = numpy.array(ordered_standard_curve)
			y_curve = numpy.array(ordered_unknown_reads)

			if self.plot is not None:
				print (x_curve)
				print (y_curve)

				self.plot.plot(x_curve, y_curve)
			
			else:
				print ("Plot needs to exist with standards before we add unknowns.")


	"""
	Shows the plot if it exists
	"""
	def show_plot(self):
		self.plot.show()

"""
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////// Methods /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
"""

"""
Selects the workbook data file holding the initial results of the experiment (raw from the imager)

	e.g. /Users/spencerdodd/Documents/Research/Broad/MPG/alamarBlue\ 
			Reading/2016_10_16\ Viral\ Titration/4.xlsx
"""
def select_data_file():
	wb = None
	while wb is None:
		print ("Please select a valid imager workbook file for analysis")
		
		try:
			root = Tkinter.Tk()
			root.withdraw()
			file_path = tkFileDialog.askopenfilename()
			
			wb = xlrd.open_workbook(file_path)
		except:
			print ("Please select a valid workbook")

	print ("Analyzing: {}".format(file_path))
	return file_path

"""
Takes in a vector each for the x-values and y-values of a dataset and returns the exponential
function of best fit. Function can then be used for display and back-calculation of viral
concentration given an RFU vlue of an unknown sample.

The function is returned in the form of 'numpy.poly1d(coefficients)', which takes in an x-value and
returns the appropriate y value on the curve.

Flag f_use determines what the returned function is based on the intended use of the function.

	"plotting"
	If the function is plotting, it returns a function whose functional use is to return y-values of
	a given x in relation to the curve of best fit

	"back-calc"
	If the function is back-calculation of concentration from an RFU value, it returns a function
	whose functional use is to return the x-value of a given y in relation to the curve of best
	fit.

"""
def get_exponential_fit(x_vector, y_vector, f_use="plotting"):

	"""
	The points used for displaying the plot are in the proper expected order
	with the x-values representing the standard concentrations and the y-value representing
	the fluorescence reading values.
	"""
	if f_use == "plotting":
		# calculate the exponential fit
		logx = numpy.log(x_vector)
		logy = numpy.log(y_vector)
		coeffs = numpy.polyfit(logx, logy, deg=3)
		exp = numpy.poly1d(coeffs)

		return exp

	"""
	The points that are used to calculate the concentrations of unknown have the x-values as 
	the fluorescence reading and the y-values as the standard concentrations (reversed).

	This is due to the fact that our return function is in a form that returns the y-value of
	a given x. Functionally, we want to get the x-value of a functional y. By reversing our x
	and y arrays that are used for fitting, we achieve this goal.
	"""
	if f_use == "back-calc":
		# reverse the x and y vectors for equation calculation
		calc_x_vector = y_vector
		calc_y_vector = x_vector

		# calculate the exponential fit
		logx = numpy.log(calc_x_vector)
		logy = numpy.log(calc_y_vector)
		coeffs = numpy.polyfit(logx, logy, deg=3)
		exp = numpy.poly1d(coeffs)

		return exp

	else:
		print ("Please enter a valid flag for the use of the returned function")

def main():
	#data_file = select_data_file()
	data_file = "/Users/spencerdodd/Documents/Research/Broad/MPG/alamarBlue Reading/2016_10_16 Viral Titration/6.xlsx"
	plate = Plate(data_file, 1900000000)
	plate.plot_standards()
	plate.plot_unknowns_with_standard_curve()
	plate.show_plot()

	"""
	data_files = [
		"/Users/spencerdodd/Documents/Research/Broad/MPG/alamarBlue Reading/2016_10_16 Viral Titration/4.xlsx",
		"/Users/spencerdodd/Documents/Research/Broad/MPG/alamarBlue Reading/2016_10_16 Viral Titration/5.xlsx",
		"/Users/spencerdodd/Documents/Research/Broad/MPG/alamarBlue Reading/2016_10_16 Viral Titration/6.xlsx",
	]
	

	for data_file in data_files:
		plate = Plate(data_file, 1900000000)
		plate.plot_standards()
		#plate.plot_unknowns_on_standard_curve()
		plate.plot_unknowns_with_standard_curve()
		plate.show_plot()
	"""

if __name__ == "__main__":
	main()




