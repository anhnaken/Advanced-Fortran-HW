DESCRIPTION OF THIS PROGRAM:

This program is utilizing the Crank-Nicolson method to solve the time-dependent Schrodinger equation. 

WHAT IS ATTACHED TO THIS FILE? And also things to look for (by the user) and some of the author's note

1. "widthsvstime.pdf"
This is a pdf file that has the graph of the spread that is evolved over time. The 
data input is written in the title but the USER can change the values that are written in the program by changing the "input6.dat"

2. "expectationvalues.pdf" 
This is a pdf file that has a graph of the expectation value evolved over time.
Once again the data input is written on the title but the USER can experiment and try different values to their best interest. 

NOTE: It can be tricky to get good values, I've noticed that even messing with values by a little bit can change the graph drastically... An example would be changing the width value. But then that could be fixed with playing around with the other values, basically the USER have to be the one to decide if it works. 

3. How to run this program: 
Simply write:
		gfortran -o Luu_proj6.x Luu_proj6.f90 matinv.f 
. The executable does not have to end in ".x" it is a matter of taste
. Notice that we have to include the "matinv.f" because there is a matrix inversion in this problem
. attached to this file is the author's own executable file: "Luu_proj6.x


4. Changing input files alongside with the executable file"

After the USER finishes putting values into "input6.dat" 
simply write: "executable file" <input6.dat to skip time and doing user input (much more efficient.

5. A lot of data 

In this file USER will find a lot of files that have .dat extensions and the author would like to take some time to clear up any confusion if there is any. 
. "analytic.dat" : This is the data file that holds the numerical results for the analytic answer for the width of the given gaussian wave function 

. "expectation.dat" : this is the data file that holds the results for the expectation values that was later put on to graph

."width.dat": this is the data file that holds the results for the calculated widths 

."xvsProbx.dat" this is the data file that writes the values of x and norm of PSI, hence the name x vs probability. 

."normalization.dat" this is the data file that holds the results for the normalization values. 






