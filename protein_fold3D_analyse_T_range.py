import random
import math
import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib import interactive
from matplotlib.ticker import NullFormatter
from mpl_toolkits import mplot3d

interactive(True)

def generate_protein(protein_length, monomer_number):
   # generate 1 dimensional array of length "protein_length"
   # each element of the array is selected randomly from 1:monomer_number
   A = np.random.randint(1, monomer_number+1, [1, protein_length]);
   
   init_x = [1] * protein_length;
   init_y = [1] * protein_length;
   init_z = [1] * protein_length;
   
   init_x[1-1] = 5;
   init_y[1-1] = 5;
   init_z[1-1] = 5;
   
   for m in range(2,protein_length+1):
        case = np.random.randint(1,4);
        if case == 1:
            init_x[m-1] = init_x[m-1-1] + 1;
            init_y[m-1] = init_y[m-1-1];
            init_z[m-1] = init_z[m-1-1];
        elif case == 2:
            init_x[m-1] = init_x[m-1-1];
            init_y[m-1] = init_y[m-1-1] + 1;
            init_z[m-1] = init_z[m-1-1];
        elif case ==3:
            init_x[m-1] = init_x[m-1-1];
            init_y[m-1] = init_y[m-1-1];
            init_z[m-1] = init_z[m-1-1] + 1;
 # concatenate A with x and y coordinates to generate a protein
 # represented as an array of monomers with associated coordinates
   protein = np.vstack((A,init_x,init_y,init_z))
   print("Protein are: ",protein)

   return protein;
   
def fold_protein(protein, T,J, number_of_runs):
    
    E_of_protein = [0]*number_of_runs;
    L_of_protein = [0]*number_of_runs;
    S_of_protein = [0]*number_of_runs;
    
    protein_length = len(protein[0]);
    E_current = protein_energy(protein, J, protein_length);
    # Choose a link at random and see if it can be moved
    for step in range(number_of_runs):
            print("Fold step: ",step)
            # Choose legal move
            copy_protein = propose_move(protein);
            E_after_move = protein_energy(copy_protein, J, protein_length);
            E_current = protein_energy(protein, J, protein_length);
            delta_E = E_after_move - E_current;
            
            if delta_E < 0: # If energy decreases, always move
              protein = copy_protein;
              E_current = E_after_move;
            else:
                boltzmann_factor = math.exp(-delta_E / T);
                if boltzmann_factor > np.random.rand():
                  protein = copy_protein;
                  E_current = E_after_move;
                  
            E_of_protein[step] = E_current;
            S_of_protein[step] = protein_size(protein,protein_length);
            L_of_protein[step] = length_end_to_end(protein, protein_length);
        
    return [E_of_protein, S_of_protein, L_of_protein, protein]

def length_end_to_end(protein, protein_length):
    # This function measures the distance between the first and last
    # monomers on the protein
    x_2 = protein[2-1][1-1];
    y_2 = protein[3-1][1-1];
    x_1 = protein[2-1][protein_length-1];
    y_1 = protein[3-1][protein_length-1];
    length = math.sqrt(((x_2 - x_1)**2) + ((y_2 - y_1)**2));  
    return length 
    
def propose_move( protein ):
    
    x_new = 0;
    y_new = 0;
    x_new2 = 0;
    y_new2 = 0;
    
    link_number = np.random.randint(1,protein_length+1);   # pick random monomer on chain
    if link_number == 1 | link_number == protein_length:  # Monomers at end of chain
        if link_number == 1:
            n = 2;
        else:
            n = protein_length - 1;
        direction = math.ceil(np.random.rand()*6);
        [x_new,y_new,z_new] = find_direction(direction,protein,n);
    
    elif (link_number + 2 <= protein_length) & (distance(protein, link_number + 1, link_number - 2) == 1):
         if protein[2-1][link_number - 1 - 1] != protein[2-1][link_number + 2 - 1]:
              x_new = protein[2-1][link_number - 1 - 1];
              x_new2= protein[2-1][link_number + 2 - 1];
              if protein[3-1][link_number-1] == protein[3-1][link_number - 1 - 1]:
                  z_new = protein[4-1][link_number - 1 - 1];
                  z_new2= protein[4-1][link_number + 2 - 1];
                  if np.random.rand() > 0.5:
                      y_new = protein[3-1][link_number - 1-1] + 1;
                      y_new2= protein[3-1][link_number + 2-1] + 1;
                  else:
                      y_new = protein[3-1][link_number - 1-1] - 1;
                      y_new2= protein[3-1][link_number + 2-1] - 1;
              elif protein[4-1][link_number-1] == protein[4-1][link_number - 1-1]:
                   y_new = protein[3-1][link_number - 1-1];
                   y_new2= protein[3-1][link_number + 2-1];
                   
                   if np.random.rand() > 0.5:
                        z_new = protein[4-1][link_number - 1-1] + 1;
                        z_new2= protein[4-1][link_number + 2-1] + 1;
                   else:
                        z_new = protein[4-1][link_number - 1-1] - 1;
                        z_new2= protein[4-1][link_number + 2-1] - 1;
         elif protein[3-1][link_number - 1-1] != protein[3-1][link_number + 2-1]:
             y_new = protein[3-1][link_number - 1-1];
             y_new2= protein[3-1][link_number + 2-1];
             
             if protein[2-1][link_number-1] == protein[2-1][link_number - 1-1]:
                  z_new = protein[4-1][link_number - 1-1];
                  z_new2= protein[4-1][link_number + 2-1];
                  
                  if np.random.rand() > 0.5:
                       x_new = protein[2-1][link_number - 1-1] + 1;
                       x_new2= protein[2-1][link_number + 2-1] + 1;
                  else:
                       x_new = protein[2-1][link_number - 1-1] - 1;
                       x_new2= protein[2-1][link_number + 2-1] - 1;
             elif protein[4-1][link_number-1] == protein[4-1][link_number - 1-1]:
                 x_new = protein[2-1][link_number - 1-1];
                 x_new2= protein[2-1][link_number + 2-1];
                 
                 if np.random.rand() > 0.5:
                    z_new = protein[4-1][link_number - 1-1] + 1;
                    z_new2= protein[4-1][link_number + 2-1] + 1;
                 else:
                    z_new = protein[4-1][link_number - 1-1] - 1;
                    z_new2= protein[4-1][link_number + 2-1] - 1;
             elif protein[4-1][link_number - 1-1] != protein[4-1][link_number + 2-1]:
                 z_new = protein[4-1][link_number - 1-1];
                 z_new2= protein[4-1][link_number + 2-1];
                 
                 if protein[2-1][link_number-1] == protein[2-1][link_number - 1-1]:
                     y_new = protein[3-1][link_number - 1-1];
                     y_new2= protein[3-1][link_number + 2-1];
                     
                     if np.random.rand() > 0.5:
                         x_new = protein[2-1][link_number - 1-1] + 1;
                         x_new2= protein[2-1][link_number + 2-1] + 1;
                     else:
                         x_new = protein[2-1][link_number - 1-1] - 1;
                         x_new2= protein[2-1][link_number + 2-1] - 1;
                         
                 elif protein[3-1][link_number-1] == protein[3-1][link_number - 1-1]:
                      x_new = protein[2-1][link_number - 1-1];
                      x_new2= protein[2-1][link_number + 1-1];
                      
                      if np.random.rand() > 0.5:
                          y_new = protein[3-1][link_number - 1-1] + 1;
                          y_new2= protein[3-1][link_number + 2-1] + 1;
                      else:
                          y_new = protein[3-1][link_number - 1-1] - 1;
                          y_new2= protein[3-1][link_number + 2-1] - 1;
    elif (link_number - 2 > 0) & (distance(protein, link_number + 1, link_number - 2) == 1) : 
        if protein[2-1][link_number + 1-1] != protein[2-1][link_number - 2-1]:
            x_new = protein[2-1][link_number + 1-1];
            x_new2= protein[2-1][link_number - 2-1];   
            
            if protein[3-1][link_number-1] == protein[3-1][link_number + 1-1]:
                z_new = protein[4-1][link_number + 1-1];
                z_new2= protein[4-1][link_number - 2-1];   
                
                if np.random.rand() > 0.5:
                    y_new = protein[3-1][link_number + 1-1] + 1;
                    y_new2= protein[3-1][link_number - 2-1] + 1;
                else:
                    y_new = protein[3-1][link_number + 1-1] - 1;
                    y_new2= protein[3-1][link_number - 2-1] - 1;
                    
            elif protein[4-1][link_number-1] == protein[4-1][link_number + 1-1]:
                y_new = protein[3-1][link_number + 1-1];
                y_new2= protein[3-1][link_number - 2-1];
                
                if np.random.rand() > 0.5:
                    z_new = protein[4-1][link_number + 1-1] + 1;
                    z_new2= protein[4-1][link_number - 2-1] + 1;
                else:
                    z_new = protein[4-1][link_number + 1-1] - 1;
                    z_new2= protein[4-1][link_number - 2-1] - 1; 
                    
        elif protein[3-1][link_number + 1-1] != protein[3-1][link_number - 2-1]:
            y_new = protein[3-1][link_number + 1-1];
            y_new2= protein[3-1][link_number - 2-1];
            
            if protein[2-1][link_number-1] == protein[2-1][link_number + 1-1]:
                z_new = protein[4-1][link_number + 1-1];
                z_new2= protein[4-1][link_number - 2-1];
                
                if np.random.rand() > 0.5:
                    x_new = protein[2-1][link_number + 1-1] + 1;
                    x_new2= protein[2-1][link_number - 2-1] + 1;
                else:
                    x_new = protein[2-1][link_number + 1-1] - 1;
                    x_new2= protein[2-1][link_number - 2-1] - 1;
                    
            elif protein[4-1][link_number-1] == protein[4-1][link_number + 1-1]:
                x_new = protein[2-1][link_number + 1-1];
                x_new2= protein[2-1][link_number - 2-1];
                
                if np.random.rand() > 0.5:
                    z_new = protein[4-1][link_number + 1-1] + 1;
                    z_new2= protein[4-1][link_number - 2-1] + 1;
                else:
                    z_new = protein[4-1][link_number + 1-1] - 1;
                    z_new2= protein[4-1][link_number - 2-1] - 1;
                    
        elif protein[4-1][link_number + 1-1] != protein[4-1][link_number - 2-1]:
             z_new = protein[4-1][link_number + 1-1];
             z_new2= protein[4-1][link_number - 2-1];
             
             if protein[2-1][link_number-1] == protein[2-1][link_number + 1-1]:
                y_new = protein[3-1-1][link_number + 1-1];
                y_new2= protein[3-1][link_number - 2-1];
                
                if np.random.rand() > 0.5:
                    x_new = protein[2-1][link_number + 1-1] + 1;
                    x_new2= protein[2-1][link_number - 2-1] + 1;
                else:
                    x_new = protein[2-1][link_number + 1-1] - 1;
                    x_new2= protein[2-1][link_number - 2-1] - 1;
                    
             elif protein[3-1][link_number-1] == protein[3-1][link_number + 1-1]:
                 x_new = protein[2-1][link_number + 1-1];
                 x_new2= protein[2-1][link_number - 2-1];
                 
                 if np.random.rand() > 0.5:
                    y_new = protein[3-1][link_number + 1-1] + 1;
                    y_new2= protein[3-1][link_number - 2-1] + 1;
                 else:
                    y_new = protein[3-1][link_number + 1-1] - 1;
                    y_new2= protein[3-1][link_number - 2-1] - 1;
                    
        check1 = not(site_occupied(x_new, y_new, z_new, protein));
        check2 = not(site_occupied(x_new2, y_new2, z_new2, protein));
        
        if  check1 & check2:
            protein[2-1][link_number-1] = x_new;
            protein[3-1][link_number-1] = y_new;
            protein[4-1][link_number-1] = z_new;
            protein[2-1][link_number - 1-1] = x_new2;
            protein[3-1][link_number - 1-1] = y_new2;
            protein[4-1][link_number - 1-1] = z_new2;    
    
    else:
        n = link_number;
        x_new = protein[2-1][n-1];
        y_new = protein[3-1][n-1];
        z_new = protein[4-1][n-1];   
        
        case = math.ceil(np.random.rand()*8);
        if case == 1:
            if ((protein[2-1][n-1]) == protein[2-1][n + 1 - 1]) & (protein[3-1][n - 1] == protein[3-1][n - 1 - 1]):
                x_new = protein[2-1][n - 1 -1];
                y_new = protein[3-1][n + 1-1];
                z_new = protein[4-1][n-1];
 
        if case == 2:
            if ((protein[2-1][n-1]) == protein[2-1][n - 1-1]) & (protein[3-1][n-1] == protein[3-1][n + 1-1]):
                x_new = protein[2-1][n + 1 -1];
                y_new = protein[3-1][n - 1-1];
                z_new = protein[4-1][n-1];
 
        if case == 3:
            if ((protein[4-1][n-1]) == protein[4-1][n + 1-1]) & (protein[2-1][n-1] == protein[2-1][n - 1-1]):
                x_new = protein[2-1][n + 1-1];
                y_new = protein[3-1][n-1];
                z_new = protein[4-1][n - 1-1];
 
        if case == 4:
            if ((protein[4-1][n-1]) == protein[4-1][n - 1-1]) & (protein[2-1][n-1] == protein[2-1][n + 1-1]):
                x_new = protein[2-1][n - 1-1];
                y_new = protein[3-1][n-1];
                z_new = protein[4-1][n + 1-1];
 
        if case == 5:
            if ((protein[4-1][n-1]) == protein[4-1][n + 1-1]) & (protein[3-1][n-1] == protein[3-1][n - 1-1]):
                x_new = protein[2-1][n-1];
                y_new = protein[3-1][n + 1-1];
                z_new = protein[4-1][n - 1-1];
 
        if case == 6:
            if ((protein[4-1][n-1]) == protein[4-1][n - 1-1]) & (protein[3-1][n-1] == protein[3-1][n + 1-1]):
                x_new = protein[2-1][n-1];
                y_new = protein[3-1][n - 1-1];
                z_new = protein[4-1][n + 1-1];
    print(x_new)        
       # [x_new,y_new,z_new] = find_protein(case,protein,n,x_new,y_new,z_new);
                 
    check3 = not(check_stretch(protein, protein_length, link_number, x_new, y_new, z_new));
    check4 = not(site_occupied(x_new, y_new, z_new, protein));
    if  check3 & check4:
        protein[2-1][link_number-1] = x_new;
        protein[3-1][link_number-1] = y_new;
        protein[4-1][link_number-1] = z_new;
          
    return protein
    
def distance(protein,n1,n2):
    distance_value = math.sqrt(((protein[2-1][n1-1] - protein[2-1][n2-1])**2) + ((protein[3-1][n1-1] - protein[3-1][n2-1])**2) + ((protein[4-1][n1-1] - protein[4-1][n2-1])**2));
    return distance_value 
    
def site_occupied(x, y,z, protein):
    
    match_x = [i for i,j in enumerate(protein[2-1]) if j == x]
    match_y = [i for i,j in enumerate(protein[3-1]) if j == y]
    match_z = [i for i,j in enumerate(protein[4-1]) if j == z]
    
    match_xy = [val for val in match_x if val in match_y];
    match_xyz = [val for val in match_xy if val in match_z];
    if not match_xyz:
        match_out = False
    else:
        match_out = True

    return match_out
    
def find_direction(direction,protein,n):
    x_new = 0;
    y_new = 0;
    z_new = 0;
    if direction == 1:
        x_new = protein[2-1][n-1] + 1;
        y_new = protein[3-1][n-1];
        z_new = protein[4-1][n-1];
 
    elif direction == 2:
        x_new = protein[2-1][n-1] - 1;
        y_new = protein[3-1][n-1];
        z_new = protein[4-1][n-1];
 
    elif direction == 3:
        x_new = protein[2-1][n-1];
        y_new = protein[3-1][n-1] + 1;
        z_new = protein[4-1][n-1];
 
    elif direction == 4:
        x_new = protein[2-1][n-1];
        y_new = protein[3-1][n-1] - 1;
        z_new = protein[4-1][n-1];
 
    elif direction == 5:
        x_new = protein[2-1][n-1];
        y_new = protein[3-1][n-1];
        z_new = protein[4-1][n-1] + 1;
 
    elif direction == 6:
        x_new = protein[2-1][n-1];
        y_new = protein[3-1][n-1];
        z_new = protein[4-1][n-1] - 1;
        
    # Execute the function
    return [x_new,y_new,z_new]
 
def check_stretch(protein, protein_length, link_number, x_new, y_new, z_new):
    
    stretched = False;
    # Test for inner links
    if (link_number > 1) & (link_number < protein_length):
        x_left = protein[2-1][link_number-1-1];
        y_left = protein[3-1][link_number-1-1];
        z_left = protein[4-1][link_number-1-1];
        x_right = protein[2-1][link_number+1-1];
        y_right = protein[3-1][link_number+1-1];
        z_right = protein[4-1][link_number+1-1];
        
        # record distances from neighbouring monomers 
        x_from_left = abs(x_new - x_left);
        x_from_right= abs(x_new - x_right);
        y_from_left = abs(y_new - y_left);
        y_from_right= abs(y_new - y_right);
        z_from_left = abs(z_new - z_left);
        z_from_right= abs(z_new - z_right);
        
        # If distance chosen location for monomer to move to is greater than 1
        # then there is a stretch
    
        if (x_from_left + y_from_left + z_from_left > 1) | (x_from_right + y_from_right + z_from_right > 1):
            stretched = True;
            
    else:
        if link_number == protein_length:    # rightmost protein
            x_nearest = protein[2-1][link_number -1-1];
            y_nearest = protein[3-1][link_number -1-1];
            z_nearest = protein[4-1][link_number -1-1];
        else:                                # leftmost protein
            x_nearest = protein[2-1][2-1];
            y_nearest = protein[3-1][2-1];
            z_nearest = protein[4-1][2-1];
    
        x_dist = abs(x_new - x_nearest);
        y_dist = abs(y_new - y_nearest);
        z_dist = abs(z_new - z_nearest);
        
        if (x_dist + y_dist + z_dist > 1):
            stretched = True;
    
    return stretched
     
def  protein_energy(protein, J, protein_length):
    
     total_energy = 0;
    
     for monomer_num in range(1,protein_length+1):
         
         x_neighbour = protein[2-1][monomer_num-1]+1;
         y_neighbour = protein[3-1][monomer_num-1];
         z_neighbour = protein[4-1][monomer_num-1];
         energy = monomer_interaction_energy(x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;        
          
        #choose neighbour below
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1]-1;
         z_neighbour = protein[4-1][monomer_num-1];
         energy = monomer_interaction_energy (x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;       
        
        # choose neighbour  left
         x_neighbour = protein[2-1][monomer_num-1]-1;
         y_neighbour = protein[3-1][monomer_num-1];
         z_neighbour = protein[4-1][monomer_num-1];
         energy = monomer_interaction_energy (x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;        
        
        # direction must be above
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1]+1;
         z_neighbour = protein[4-1][monomer_num-1];
         energy = monomer_interaction_energy (x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;
         
         # direction must be forward
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1];
         z_neighbour = protein[4-1][monomer_num-1] +1;
         energy = monomer_interaction_energy (x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;
         
         # direction must be backward
         x_neighbour = protein[2-1][monomer_num-1];
         y_neighbour = protein[3-1][monomer_num-1];
         z_neighbour = protein[4-1][monomer_num-1]-1;
         energy = monomer_interaction_energy (x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J); # This will check if
        # occupied and if so, calculate the interaction energy
         total_energy = total_energy + energy;
         
    # Since monomer interactions have been double counted, energy
    # calculated is twice as high as it should be
     return total_energy
     
def monomer_interaction_energy(x_neighbour, y_neighbour,z_neighbour, protein, monomer_num, J):
# Calculates the interaction energy between a monomer in a protein and a
# neighbouring monomer if the neighbouring monomer exists. If the input
# coordinates of a potential neighbour do not corrospond to an existing
# monomer, return 0
    energy = 0;
    if site_occupied(x_neighbour, y_neighbour,z_neighbour,protein):
        
       find_x = func_find(protein[2-1],lambda x: x == x_neighbour);
       find_y = func_find(protein[3-1],lambda y: y == y_neighbour);
       find_z = func_find(protein[4-1],lambda z: z == z_neighbour);
       
       intersect_xy = [val for val in find_x if val in find_y];
       interscet_xyz = [val for val in intersect_xy if val in find_z];
       neighbour = interscet_xyz; 
        # The interaction energy is only a,n issue for monomers which are
        # not linked on the protein chain
       if not neighbour:
           k = 0;
       else:
           k =[];
           for k in neighbour:
               print()
       u = np.subtract(k+1,monomer_num)
       if (abs(u) > 1):
            # Use J matrix to find interaction energy between the two
            # monomers on the chain
            energy = J[int(protein[1-1][monomer_num-1])-1][int(protein[1-1][k])-1];
            
    return energy

def func_find(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]
    
def protein_size(protein,protein_length):
    
    # Initial center of mass coordinates set to 0
    c_x = 0;
    c_y = 0;
    c_z = 0;
    
    # Initial ruggedness set to 0;
    s = 0;
    
    for monomer_num in range(1,protein_length+1): 
        c_x += protein[2-1];
        c_y += protein[3-1];
        c_z += protein[4-1];
        
    np.divide(c_x,protein_length)
    np.divide(c_y,protein_length)
    np.divide(c_z,protein_length)
    
    for monomer_num in range(1,protein_length+1):
        s += (((protein[2-1] - c_x)**2) + ((protein[3-1] - c_y)**2) + ((protein[4-1] - c_z)**2));

    return (s / protein_length)


##### Main Block ##########
monomer_number = 20;
A = np.random.rand(monomer_number,monomer_number);
random.shuffle(A)
A = np.multiply(A,2)
J = A - 4;

T_range = np.arange(9, 10); # temperature in kelvin
protein_length = 15;
number_of_runs = 5000;

init_protein = generate_protein(protein_length,monomer_number);

E_temp = []
L_temp = []
S_temp = []
E_vs_T = [0] * len(T_range);
L_vs_T= [0] * len(T_range);
S_vs_T= [0] * len(T_range);

step = 0
for T in range(1,len(T_range)+1):
    print('Main Step:', step)
    [E_temp, S_temp, L_temp, final_protein] = fold_protein(init_protein, T,J, number_of_runs);
    E_vs_T[step] = np.mean(E_temp);
    L_vs_T[step] = np.mean(L_temp);
    S_vs_T[step] = np.mean(S_temp);
    step = step + 1;

# plot graphs

plt.figure(1)

# display energy of protein at each step
plt.subplot(221)
plt.plot(E_vs_T)
plt.xlabel('Monte Carlo steps');
plt.ylabel('Energy')
plt.legend ('Energy vs time');  

# display "end to end" length of protein
plt.subplot(222)
plt.plot(L_vs_T)
plt.xlabel('Monte Carlo steps');
plt.ylabel('Length')
plt.legend('End to end length')

# display "size of protein" 
plt.subplot(223)
plt.plot(S_vs_T)
plt.xlabel('Monte Carlo steps');
plt.ylabel('Protein Size')
plt.legend('size')

plt.figure(2)
# display protein lattice
ax = plt.axes(projection='3d')
# Data for a three-dimensional line

x = final_protein[2-1];
y = final_protein[3-1];
z = final_protein[4-1];
ax.plot3D(x, y, z, 'red')
#ax.set_xlim(0, 30); ax.set_ylim(0, 30)
plt.legend('Protein Lattice')

plt.gca().yaxis.set_minor_formatter(NullFormatter())
# Adjust the subplot layout, because the logit one may take more space
# than usual, due to y-tick labels like "1 - 10^{-3}"
plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.25,wspace=0.35)

plt.show()