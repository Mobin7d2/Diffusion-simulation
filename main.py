import numpy as np
import matplotlib.pyplot as plt
from random import randint

# General variables
# Data corresponding to each molecule location is stored in the following dictionaries
Gas_1 = {}
Gas_2 = {}


def Random_walk_update_1(Data_1, Data_2, boundary):
    #This function updates the location of each molecule of Gas 1.
    #It should be noted that data contains coordinate of each molecule as [x, y] .
    #Gas 1 can freely move onto every empty vertex
    #boundary is of the form [x, y] where x and y are the considered as the boundary of the container
    Molecule = list(Data_1.keys())
    for i in Molecule:
        # 1 right 2 left 3 up 4 down
        Next_step = randint(1, 4)
        if Next_step == 1 :
            # One step to right
            Data_1[i][0] += 1
        elif Next_step == 2:
            # One step to left
            Data_1[i][0] -= 1
        elif Next_step == 3:
            # One step up
            Data_1[i][1] += 1
        else:
            # One step down
            Data_1[i][1] -= 1


        ### No more than one
        ## In this section we check if more than one molecule occupy a single vertex or not
        # This variable contians position of all Molecules
        Values = list(Data_1.values()) + list(Data_2.values())
        counter = 0
        for position in Values:
            if position == Data_1[i]:
                counter += 1

            if counter > 1 :
                #This section practically redo the same calculation for this molecule
                i -= 1
                break
        

        ### Boundary check:
        # This section checks if the molecule has exited the container or not
        if not counter > 1:
            if Data_1[i][0] == boundary[0]:
                Data_1[i][0] -= 1
                if i == 0:
                    Random_walk_update_1(Data_1, Data_2, boundary)
                else: 
                    i -= 1
                
            elif Data_1[i][0] == -boundary[0]:
                Data_1[i][0] += 1
                if i == 0:
                    Random_walk_update_1(Data_1, Data_2, boundary)
                    
                else: 
                    i -= 1
                

            if Data_1[i][1] == boundary[1]:
                Data_1[i][1] -= 1
                if i == 0:
                    Random_walk_update_1(Data_1, Data_2, boundary)
                    
                else: 
                    i -= 1
                
            elif Data_1[i][1] == -boundary[1]:
                Data_1[i][1] += 1
                if i == 0:
                    Random_walk_update_1(Data_1, Data_2, boundary)
                else: 
                    i -= 1
                

        

            

def Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1):
    #This function updates the location of each molecule.
    #It should be noted that data contains coordinate of each molecule as [x, y] .
    Molecule = list(Data_2.keys())
    for i in Molecule:
        # 1 right 2 left 3 up 4 down
        Next_step = randint(1, 4)
        if Next_step == 1 :
            ### boundary check
            if Data_2[i][0] + 1 == - boundary_container_1[0] or Data_2[i][0] + 1 >= -boundary_container_1[0]:
                i -= 1
                break
            # One step to right
            Data_2[i][0] += 1

            
        elif Next_step == 2:
            ### boundary check
            if Data_2[i][0] - 1 == boundary_container_1[0] or Data_2[i][0] - 1 >= boundary_container_1[0]:
                i -= 1
                break

            # One step to left
            Data_2[i][0] -= 1
        elif Next_step == 3:
            ### boundary check
            if Data_2[i][1] + 1 ==  -boundary_container_1[1]or Data_2[i][1] + 1 >= -boundary_container_1[1]:
                i -= 1
                break

            # One step up
            Data_2[i][1] += 1
        else:
            ### boundary check
            if Data_2[i][1] - 1 ==  boundary_container_1[1] or Data_2[i][1] - 1 >= boundary_container_1[1]:
                i -= 1
                break

            # One step down
            Data_2[i][1] -= 1

        ### No more than one 
        ## In this section we check if more than one molecule occupy a single vertex or not
        Values = list(Data_2.values()) + list(Data_1.values())
        counter = 0
        for position in Values:
            if position == Data_2[i]:
                counter += 1

            if counter > 1 :
                #This section practically redo the same calculation for this molecule
                i -= 1
                break
        


        ### Boundary check:
        # This section checks if the molecule has exited the container or not
        if not counter > 1:
            # This parameter ends this process the moment one of the conditions is violated
            flag = False
            # X component check
            if Data_2[i][0] == boundary[0]:
                Data_2[i][0] -= 1
                if i == 0:
                    Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                    
                else: 
                    i -= 1
                    flag = True
            elif Data_2[i][0] == -boundary[0]:
                Data_2[i][0] += 1
                if i == 0:
                    Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                    
                else: 
                    i -= 1
                    flag = True

            # Y component check
            if Data_2[i][1] == boundary[1]:
                Data_2[i][1] -= 1
                if i == 0:
                    Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                    
                else: 
                    i -= 1
                    flag = True
            elif Data_2[i][1] == -boundary[1]:
                Data_2[i][1] += 1
                if i == 0:
                    Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                    
                else: 
                    i -= 1
                    flag = True


            # Is it inside container 1 ?
            if not flag:
                # X component check
                if Data_2[i][0] == boundary_container_1[0]:
                    Data_2[i][0] += 1
                    if i == 0:
                        Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                        
                    else: 
                        i -= 1
                elif Data_2[i][0] == -boundary_container_1[0]:
                    Data_2[i][0] -= 1
                    if i == 0:
                        Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                        
                    else: 
                        i -= 1
                # Y componenet check
                if Data_2[i][1] == boundary_container_1[1]:
                    Data_2[i][1] += 1
                    if i == 0:
                        Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                        
                    else: 
                        i -= 1
                        flag = True
                elif Data_2[i][1] == -boundary_container_1[1]:
                    Data_2[i][1] -= 1
                    if i == 0:
                        Random_walk_update_2(Data_1, Data_2, boundary, boundary_container_1)
                        
                    else: 
                        i -= 1
                        flag = True

        




def Initial_coordinate(width, height, N):
    # This function returns initial coordinates of N molecule
    Coordinates = []
    # Key values of this variable is their x and each value is y of the vertex regarding
    #  that x which has not yet been taken
    Allowed_coordinates = {}

    # width is given in form of [x1, x2]. Width of this rigion is considered
    # as x1 < x < x2. The same logic is used for heigth.
    if width[0] != 0:
        # Molecule 2
        Allowed_x = list(range(-width[1], width[1]+1))
        Allowed_y = list(range(-height[1], -height[0])) + list(range(height[0]+1, height[1]+1))
        All_y = list(range(-height[1], height[1]+1))
    else:
        Allowed_x = list(range(-width[1], -width[0]+1)) + list(range(width[0]+1, width[1]+1)) 
        Allowed_y = list(range(-height[1], -height[0]+1)) + list(range(height[0]+1, height[1]+1))
    
    
    for j in Allowed_x:
        # Setting up the new variable Allowed_coordinates
        if -width[0] <= j <= width[0]:
            Allowed_coordinates[j] = Allowed_y
        else:
            Allowed_coordinates[j] = list(range(-height[1], height[1]+1))

    

    # Note that N has to be smaller than the number of dots in our grid
    for i in range(N):
        Allowed_coordinates_x = list(Allowed_coordinates.keys())
        X_random_choice_index = randint(0, len(Allowed_coordinates_x)-1)
        X_random_choice = Allowed_coordinates_x[X_random_choice_index]

        # possible choices for the x taken
        Y_possible_choice = Allowed_coordinates[X_random_choice]
        Y_random_choice_index = randint(0, len(Y_possible_choice)-1)
        Coordinates.append([Allowed_coordinates_x[X_random_choice_index],
                             Y_possible_choice[Y_random_choice_index]])


        # Updating the possibilities so that no two molecule would occupy the same place.
        Y_new_possible_choice = []
        for j in Y_possible_choice:
            if j != Y_possible_choice[Y_random_choice_index]:
                Y_new_possible_choice.append(j)

        Allowed_coordinates.update({X_random_choice:Y_new_possible_choice})
        if Allowed_coordinates[X_random_choice] == []:
            # Deletes empty columns of x
            Allowed_coordinates.pop(X_random_choice)
            
    
    
    return Coordinates


def X_Y_molecule(Gas):
    #This function converts data of the dictionary containing coordinates of each molecules to lists!
    data = list(Gas.values())
    Coordinates = [[], []]
    for element in data:
        Coordinates[0].append(element[0])
        Coordinates[1].append(element[1])

    return Coordinates


def Entropy(Molecules_position, Square_width, N_walkers):
    #We need to make a grid on the square b
    Number_cells = Square_width//2 
    Cells = np.zeros((Number_cells, Number_cells))
    
    Positions = list(Molecules_position.values())
    for pos in Positions:
        x = pos[0] // Number_cells
        y = pos[1] // Number_cells
        Cells[x, y] += 1
        
    probability = Cells/N_walkers
    probability = probability[probability>0] #Excluding zeros

    k_b = 1.38*10**(-23)
    entropy = -k_b * np.sum(probability * np.log(probability))
    return entropy




def main(Number_of_iteration, N1, N2):
    # a and b square (almost half length)
    a_width = 10
    a_hight = 10
    b_width = 15
    b_hight = 15
    # Initial position for gas 1
    Coordinates = Initial_coordinate([0, a_width-1], [0, a_hight-1], N1)
    for i in range(N1):
        Gas_1[i] = Coordinates[i]
    
    
    # Initial position for gas 2
    Coordinates = Initial_coordinate([a_width, b_width-1], [a_hight, b_hight-1], N2)
    for i in range(N2):
        Gas_2[i] = Coordinates[i]

    # The following section determines how the coordinate of each molecule changes over a 
    # number of iteration.
    
    Time = [i for i in range(Number_of_iterations)]
    Entropy_value = []

    #Plotting 
    
    for i in range(Number_of_iterations):
        Entropy_value.append(Entropy(Gas_1, b_width, N1))
        Random_walk_update_1(Gas_1, Gas_2, [b_width, b_hight])
        Random_walk_update_2(Gas_1, Gas_2, [b_width, b_hight], [a_width, a_hight])
        
        if i % 400 == 0 or i == 1999:
            #Plotting position of molecules     
            # Setting geometry of the graph
            f = plt.figure()
            f.set_figwidth(8)
            f.set_figheight(8)
    
            #Plotting boarders of the containers
            # plotting a line plot after changing it's width and height
            plt.plot([a_width, -a_width, -a_width, a_width, a_width], [a_hight, a_hight, -a_hight, -a_hight, a_hight], color = "black")
            plt.plot([b_width, -b_width, -b_width, b_width, b_width], [b_hight, b_hight, -b_hight, -b_hight, b_hight], color = "black")
            plt.title("Molecules' position")
            plt.scatter(X_Y_molecule(Gas_1)[0], X_Y_molecule(Gas_1)[1], color = "blue", label="Gas 1")
            plt.scatter(X_Y_molecule(Gas_2)[0], X_Y_molecule(Gas_2)[1], color = "red", label="Gas 2")
            plt.legend()
            plt.savefig("./img/"+str(i))
            #plt.show()

    Averaged_entropy_value = []
    sum = 0
    for j in range(len(Entropy_value)):
        sum+= Entropy_value[i]
        if (j+1) % 2 ==0:
            Averaged_entropy_value.append(sum/10)
            sum = 0 
    

    
    
    #Plotting entropy vs time
    f = plt.figure()
    f.set_figwidth(8)
    f.set_figheight(8)
    plt.title("Entropy Vs time")
    plt.xlabel("Time(s)")
    plt.ylabel("Entropy(J/K)")
    plt.scatter(Time, Entropy_value, color="orange")
    plt.savefig("./img/EntropyVSTime_"+str(N1)+"_"+str(N2))
    #plt.show()


Number_of_iterations = 2000
main(Number_of_iterations, 80, 0)
