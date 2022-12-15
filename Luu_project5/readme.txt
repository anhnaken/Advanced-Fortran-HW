DESCRIPTION OF THIS PROGRAM: 

This program is utilizing the method of the Runge-Kutta algorithm to solve differential equations.

WHAT IS ATTACHED TO THIS FILE? Things to Look For and maybe some of the Author's Notes:

1. "MoonandPlanetOrbit.pdf" 
This is a pdf file that has the graph of a moon orbiting a planet, 
the data input is written in the title but USER can change values with 
"mooninputs.dat" file.

2. Two planet system orbits:
USER can see different graphs of a two planet system in orbit.
USER CAN FIND THESE GRAPHS IN THE FOLDER LABELED: "graphs" 
LOOK FOR THESE GRAPHS: 
	. "TwoPlanets1.pdf"
	. "TwoPlanets2.pdf"
	. "TwoPlanets3.pdf"
	. "TwoPlanets4.pdf"

NOTE: EACH of these graphs shows a two planet system in orbit, however, they vary in input values(WHICH IS FOUND AT THE TOP AT THE TITLE), also USER can change values to check with their own regards by looking at "planetinputs.dat" (WHICH AT THE MOMENT HAS THE VALUES OF "TwoPlanets4.pdf" once again user can change to their desire). 

3. How to Run this Program:
Simply write:
		 gfotran -o Luu_projec5.x Luu_projec5.f90 
		
Luu_project5.x is the executable file. 
Of course the USER does not have to put ".x" it is a matter of choice and style.  


4. How Does One Change the Input Files With the Executable File?

	. Change the values in order that is listed in the files "planetinputs.dat"
	and "mooninputs.dat" to get the wanted orbits
	. Then, with the executable (it does not have to end with ".x" but I chose it to b	be, user can choose what makes them happy) simply put: 
	"./Luu_proj5.x <mooninputs.dat" or "./Luu_proj5.x <planetinputs.dat"

5. THINGS TO LOOK OUT FOR:

	.When entering in the input files, it can be tricky from my experience because 		changing a value (even by a little bit) can drastically change the graphs. 
	Hence, user must find the appropriate values to get the graphs that is needed.
	Which is why provided is a couple of different graphs displaying different values
	but showing that they somewhat have the same consequences (graphically speaking).
	
	.For the moon orbit("MoonandPlanetOrbit.pdf") to get that, I varied my moon 
	x-values just a little higher than the planet's x values and as for the moon's
	velocity it is larger than the planet's velocity. From there I played with the 
	numbers until I got what I wanted. The USER can try different values, of course. 

	





