import random
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import interactive
from matplotlib.ticker import NullFormatter
interactive(True)

def generate_protein(protein_length, monomer_number):
   # generate 1 dimensional array of length "protein_length"
   # each element of the array is selected randomly from 1:monomer_number
   A = np.random.randint(1, monomer_number+1, [1, protein_length]);
   
   # Initial x and y coordinates for protein
   x_offset = 10;
   y_offset = 15;
   
   init_x = np.arange(x_offset, (protein_length + 9 + 1), 1);
   init_y = np.ones((1, protein_length)) * y_offset;
   
 # concatenate A with x and y coordinates to generate a protein
 # represented as an array of monomers with associated coordinates
   protein = np.vstack((A,init_x,init_y))
   #print("Protein: ",protein)

   return protein;

def length_end_to_end(protein, protein_length):
    # This function measures the distance between the first and last monomers 
    # on the protein
    x_2 = protein[2-1][1-1];
    y_2 = protein[3-1][1-1];
    x_1 = protein[2-1][protein_length-1];
    y_1 = protein[3-1][protein_length-1];
    length = math.sqrt(((x_2 - x_1)**2) + ((y_2 - y_1)**2));  
    return length 

def check_stretch(protein, protein_length, link_number, x_new, y_new):
    
    stretched = False;
    # Test for inner links
    if (link_number > 1) & (link_number < protein_length):
        x_left = protein[2-1][link_number-1-1];
        y_left = protein[3-1][link_number-1-1];
        x_right = protein[2-1][link_number+1-1];
        y_right = protein[3-1][link_number+1-1];
        
        # record distances from neighbouring monomers 
        x_from_left = abs(x_new - x_left);
        x_from_right= abs(x_new - x_right);
        y_from_left = abs(y_new - y_left);
        y_from_right= abs(y_new - y_right);
        
        # If distance chosen location for monomer to move to is greater than 1
        # then there is a stretch which results in protein breakage
        
        if (x_from_left + y_from_left > 1) | (x_from_right + y_from_right > 1):
            stretched = True;
            
    else:
        if link_number == protein_length:    # rightmost protein
            x_nearest = protein[2-1][link_number-1-1];
            y_nearest = protein[3-1][link_number-1-1];
        else:                                # leftmost protein
            x_nearest = protein[2-1][2-1];
            y_nearest = protein[3-1][2-1];
            
        x_dist = abs(x_new - x_nearest);
        y_dist = abs(y_new - y_nearest);
    
        if (x_dist + y_dist > 1):
            stretched = True;
    
    return stretched
    
def site_occupied(x, y, protein):
    # This function determines if the chosen location is occupied by another
    # monomer or not
    
    match_x = [i for i,j in enumerate(protein[2-1]) if j == x]
    match_y = [i for i,j in enumerate(protein[3-1]) if j == y]
    match_xy = [val for val in match_x if val in match_y];
    if not match_xy:
        match_out = False
    else:
        match_out = True

    return match_out
    
def find_direction(direction,protein,n):
    # This function is used to determine the monomer position after the move
    x_new = 0;
    y_new = 0;
    if direction == 1:  # RIGHT and UP
      x_new = protein[2-1][n-1] + 1;
      y_new = protein[3-1][n-1] + 1;
      
    elif direction == 2: # RIGHT and STATIC
      x_new = protein[2-1][n-1]+1;
      y_new = protein[3-1][n-1];
 
    elif direction == 3:  # RIGHT and DOWN
      x_new = protein[2-1][n-1] + 1;
      y_new = protein[3-1][n-1] - 1;
      
    elif direction == 4: # STATIC and DOWN
      x_new = protein[2-1][n-1];
      y_new = protein[3-1][n-1]-1;
 
    elif direction == 5: # LEFT and DOWN
      x_new = protein[2-1][n-1] - 1;
      y_new = protein[3-1][n-1] - 1;
      
    elif direction == 6: # LEFT and STATIC
      x_new = protein[2-1][n-1]-1;
      y_new = protein[3-1][n-1];
 
    elif direction == 7: # LEFT and UP
      x_new = protein[2-1][n-1] - 1;
      y_new = protein[3-1][n-1] + 1;
    else:               # STATIC and UP
      x_new = protein[2-1][n-1];
      y_new = protein[3-1][n-1]+1;
      
    # Execute the function
    return [x_new,y_new]
    
def  protein_energy(protein, J, protein_length):
    # Function is used to calculate total energy of the protein
     total_energy = 0;
    
     for monomer_num in range(1,protein_length+1):
         
         x_neighbour = protein[2-1][monomer_num-1]+1;
         y_neighbour = protein[3-1][monomer_num-1];
         energy = monomer_interaction_energy(x_neighbour, y_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;        
          
        #choose neighbour below
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1]-1;
         energy = monomer_interaction_energy (x_neighbour, y_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;       
         
        # choose neighbour  left
         x_neighbour = protein[2-1][monomer_num-1]-1;
         y_neighbour = protein[3-1][monomer_num-1];
         energy = monomer_interaction_energy (x_neighbour, y_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;     
        
        # direction must be above
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1]+1;
         energy = monomer_interaction_energy (x_neighbour, y_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy  
         total_energy=total_energy+energy; 

     return total_energy
     
def monomer_interaction_energy(x_neighbour, y_neighbour, protein, monomer_num, J):
# Function is used to calculate energy at each amino acid. 
    energy = 0;
    if site_occupied(x_neighbour,y_neighbour,protein):
        
       find_x = func_find(protein[2-1],lambda x: x == x_neighbour);
       find_y = func_find(protein[3-1],lambda y: y == y_neighbour);
       
       neighbour = [val for val in find_x if val in find_y];
        # The interaction energy is only an issue for monomers which are
        # not linked on the protein chain

       if not neighbour:
           k = 0;
       else:
           k =[];
           for k in neighbour:
               print()
       u = np.subtract(k+1,monomer_num) # k+1
       if (abs(u) > 1):
            # Use J matrix to find interaction energy between the two
            # monomers on the chain
            energy = J[int(protein[1-1][monomer_num-1])-1][int(protein[1-1][k])-1]; # k           
    return energy
     
def func_find(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]
    
################# MAIN BLOCK #############################

# Initialisation block
protein_length = 15;   # Set protein length(monomer number in the chain)
number_of_runs = 5000; #Monte Carlo steps
T = 10;    # temperature in kelvin

monomer_number = 20; # There are 20 monomers occuring in nature

high_interaction = -400000;
low_interaction = -200000;

E_current = 0;
Average_energy=0;

E_of_protein = [0] * number_of_runs;
L_of_protein= [0] * number_of_runs;

# J is a 20x20 matrix of randomly assigned energy values to represent the
# interaction energies betweeen monomers

J = [20]*monomer_number
for i in range(monomer_number):
    J[i] = [0] * monomer_number
for i in range(monomer_number):
    for j in range(monomer_number):
        J[i][j] = random.randint(high_interaction,low_interaction)/100000

protein = generate_protein(protein_length, monomer_number); # ok

# Choose a link at random and see if it can be moved
#protein_length=len(protein);

for step in range(number_of_runs):
    print("Fold step: ",step)
    link_number = np.random.randint(1,protein_length+1); # pick random monomer on chain
    direction = math.ceil(np.random.rand()*8);  # pick direction denoted by number from 1 to 8
    
    [x_new,y_new] = find_direction(direction,protein,link_number); # ok
    
     # check if chosen location is occupied
    occupied = site_occupied(x_new, y_new, protein); # ok
    
    # check if chosen location causes "stretch"
    stretched = check_stretch(protein, protein_length, link_number, x_new, y_new); # ok
    #sys.exit()
   # stretched = True;
    
    # if chosen location not occupied, create copy of protein with new shape
    
    check1=not(occupied);
    check2=not(stretched);
    
    if  check1 & check2:
        copy_protein = protein;
        copy_protein[2-1][link_number-1] = x_new;
        copy_protein[3-1][link_number-1] = y_new;
        # Compare energy value of new protein shape with the old shape
        
        E_after_move = protein_energy(copy_protein, J, protein_length);
        E_current = protein_energy(protein, J, protein_length);
        
        delta_E = E_after_move - E_current;
        if delta_E <= 0: # If energy decreases, always move
             protein = copy_protein;
             E_current = E_after_move;
        else:
                # if T == 0:
                #     if delta_E != 0:
                #       protein = copy_protein;
                #       E_current = E_after_move;
                # else:
             boltzmann_factor = math.exp(-delta_E / T);
             if boltzmann_factor > np.random.rand():
                 protein = copy_protein;
                 E_current = E_after_move;
    
    E_of_protein[step] = E_current;
    L_of_protein[step] = length_end_to_end(protein, protein_length);
    Average_energy = Average_energy + E_current
    
# plot graphs
Eavg = Average_energy/number_of_runs;
print('Average_Energy: ', Eavg)

plt.figure(1)

# display energy of protein at each step
#plt.subplot(221)

plt.plot(E_of_protein)
plt.xlabel('Monte Carlo steps');
plt.ylabel('Energy')
plt.title ('Energy vs time');  

# display "end to end" length of protein
#plt.subplot(222)
plt.figure(2)
plt.plot(L_of_protein)
plt.xlabel('Monte Carlo steps');
plt.ylabel('Length')
plt.title('End to end length')

# display protein lattice
#plt.subplot(223)
plt.figure(3)
x = protein[2-1];
y = protein[3-1];
plt.plot(x, y,'-r',markersize = 5)
plt.axis([0, 30, 0, 30])
plt.title('Protein Lattice')

#plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
#plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,wspace=0.35)

plt.show()