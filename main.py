import numpy as np
import matplotlib.pyplot as plt

# Define parameters

# Here h stands for hydrophobic amino acids, while p stands for polar amino acids
protein = "hhphhphphphphhpphpph"

size = len(protein) 
num_ants = 100 #number of ants
num_iterations = 100 # number of iterations for the ACO algorithm
alpha = 1.0
beta = 5.0
evap_rate = 0.5 # Evaporation rate
Emin = 10.0 # negative of E*min value

pheromone = np.ones((2*size+2,2*size+2,4)) # Pheromone matrix initialised to 1

def checkLetter(sign): 
    if sign == 'p':  # polar
        return False
    elif sign == 'h':  # hydrophobic
        return True
    
def replaceLetter(sign): # function used to fill the board/lattice
    if sign == 'p':  # polar
        return -1
    elif sign == 'h':  # hp
        return 1

'''
This method is used to calculate score of the folding created till now. 
Score is calculated as the number of non-consecutive adjacent Hydrophobic amino-acids
Energy of a folded protein is equal to the negative of this score. 
Hence, we want the maximal score to get minimal energy.
''' 
def score(board, length=size):
    # Horizontal check
    score = 0
    for i in range(len(board)):
        for j in range(len(board[0])-1):
            if(board[i][j]==1 and board[i][j+1]==1):
                score+=1

    # Vertical check
    for i in range(len(board)-1):
        for j in range(len(board[0])):
            if(board[i][j]==1 and board[i+1][j]==1):
                score+=1

    # Removing consecutive adjacent H amino acids count from score
    for i in range(length-1):
        if(protein[i]=='h' and protein[i+1]=='h'):
            score-=1

    return score

def calc_prob(k,dir, curx, cury, board):
    nextx = curx
    nexty = cury

    if(dir==0): #right
        nexty+=1
    if(dir==1): #up
        nextx-=1
    if(dir==2): #left
        nexty-=1
    if(dir==3): #down
        nextx+=1
    
    # Return 0 if the next point is already filled
    if(board[nextx][nexty] != 0 ):
        return 0
    
    pheromone_level = pheromone[curx][cury][dir]
    board[nextx][nexty] = replaceLetter(protein[k+1])
    # The value of eta(heuristic) is derived from the research paper
    eta = 1 + score(board,k+2)
    board[nextx][nexty] = 0 # Removing the filled point for score

    # Using the ACO probability distribution formulae
    prob = (pheromone_level ** alpha) * (eta ** beta)
    return prob

    
def select_next_node(k,curx, cury, board):
    # Initialising probability array for all directions movement
    probabilities = np.zeros(4)
    for i in range(4):
        probabilities[i] = calc_prob(k,i,curx,cury,board)

    # Handling the edge case when all possible movements are blocked, so returning -1
    if(np.sum(probabilities)==0):
        return [-1,-1]
   
    probabilities = probabilities / np.sum(probabilities)
    dir = np.random.choice([0,1,2,3], p=probabilities)
    nextx = curx
    nexty = cury

    if(dir==0):
        nexty+=1
    if(dir==1):
        nextx-=1
    if(dir==2):
        nexty-=1
    if(dir==3):
        nextx+=1
    return [nextx,nexty]

'''Function for updating pheromones after all ants have completed their paths'''
def updatePheromones(pheromone, all_pathsx, all_pathsy, all_ants_score):
    pheromone *= evap_rate
    
    for i in range(len(all_pathsx)):
        # Handling the edge case of path being wrong
        if(all_pathsx[i][0]==-1):
            continue
        for j in range(len(all_pathsx[i])-1):
            curx = all_pathsx[i][j]
            cury = all_pathsy[i][j]
            nextx = all_pathsx[i][j+1]
            nexty = all_pathsy[i][j+1]
            dir = 0
            if(nexty-cury==1):
                dir = 0
            if(nextx-curx== -1):
                dir = 1
            if(nexty-cury == -1):
                dir = 2
            if(nextx-curx == 1):
                dir = 3

            pheromone[curx][cury][dir] += all_ants_score[i]/Emin
            


def aco():
    bestpolargraphx = [] # Array to represent polar nodes' x values
    bestpolargraphy = [] # y value corresponding to x values
    besthgraphx = [] # Array to represent hydrophobic nodes' x values
    besthgraphy = []
    bestgraphx = [] # Array that has all the positions
    bestgraphy = []
    bestscore = -1.0 # Variable to handle best score among all iterations
    for i in range(num_iterations):
        all_pathsx = [] # List of paths of all the ants in this interation
        all_pathsy = []
        all_ants_score = np.zeros(num_ants) # Array of score of each ant, initialized to zero
        for j in range(num_ants):
            polargraphx = []
            polargraphy = []
            hgraphx = []
            hgraphy = []
            graphx = []
            graphy = []

            #Defining the board for each ant 
            board = np.zeros((2*size+2, 2*size+2), dtype=int)

            # Putting first amino acid in middle of the board
            positionx = int(size+1)
            positiony = int(size+1)

            graphx.append(positionx)
            graphy.append(positiony)

            board[positionx][positiony] = replaceLetter(protein[0])

            if checkLetter(protein[0]):  # hp
                hgraphx.append(positionx)
                hgraphy.append(positiony)
            else:  # polar
                polargraphx.append(positionx)
                polargraphy.append(positiony)

            flag = False # to check if current ant's path was wrong/interlinking

            # Moving the ant to complte the folding of protein
            for k in range(0,size-1):
                next = select_next_node(k,positionx,positiony,board)
                nextx = next[0]
                nexty = next[1]
                if(nextx==-1 or nexty==-1):
                    flag = True
                    break
                positionx = nextx
                positiony = nexty
                graphx.append(positionx)
                graphy.append(positiony)

                board[positionx][positiony] = replaceLetter(protein[k+1])
                if checkLetter(protein[k+1]):  # hp
                    hgraphx.append(positionx)
                    hgraphy.append(positiony)
                else:  # polar
                    polargraphx.append(positionx)
                    polargraphy.append(positiony)
            
            if(flag==True):
                all_ants_score[j] = 0
                all_pathsx.append([-1])
                all_pathsy.append([-1])
                continue

            curscore = score(board)
            #print(curscore) # checker here
            all_ants_score[j] = curscore
            all_pathsx.append(graphx)
            all_pathsy.append(graphy)

            # Update bestscore if curscore is better
            if(curscore > bestscore):
                bestscore = curscore
                bestgraphx = graphx
                bestgraphy = graphy
                besthgraphx = hgraphx
                besthgraphy = hgraphy
                bestpolargraphx = polargraphx
                bestpolargraphy = polargraphy

        updatePheromones(pheromone,all_pathsx,all_pathsy,all_ants_score)    

    return bestscore, bestgraphx, bestgraphy, besthgraphx, besthgraphy, bestpolargraphx, bestpolargraphy       

            
bestscore, bestgraphx, bestgraphy, besthgraphx, besthgraphy, bestpolargraphx, bestpolargraphy  = aco()


print("Best path: ", bestgraphx)
print("Best cost: ", bestscore)

# Plotting

plt.title("Protein Folding by ACO")
plt.plot(bestgraphy, bestgraphx, 'k')
plt.plot(bestpolargraphy, bestpolargraphx, 'co', besthgraphy, besthgraphx, 'ko', markersize=12)
plt.grid(True)
plt.axis([0, 2*size+2, 0, 2*size+2])
plt.show()