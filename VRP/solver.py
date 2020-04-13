'''Code file for vehicle routing problem created for Advanced Algorithms
Spring 2020 at Olin College. These functions solve the vehicle routing problem
using an integer programming and then a local search approach. This code has
been adapted from functions written by Alice Paul.'''

import picos as pic
from picos import RealVariable, Constant, BinaryVariable
import numpy as np
from read_files import read_file_type_A, read_file_type_C
from scipy.special import comb

# Integer programming approach
def cvrp_ip(C,q,K,Q,obj=True):
    '''
    Solves the capacitated vehicle routing problem using an integer programming
    approach.

    C: matrix of edge costs, that represent distances between each node
    q: list of demands associated with each client node
    K: number of vehicles
    Q: capacity of each vehicle
    obj: whether to set objective (ignore unless you are doing local search)
    returns:
        objective_value: value of the minimum travel cost
        x: matrix representing number of routes that use each arc
    '''
    # Add in destination node
    # It should have the same distances as source and its demand equals 0
    q2 = np.append(q, 0)
    col = C[:,0]
    C2 = np.hstack((C, np.atleast_2d(col).T))
    row = C2[0,:]
    C3 = np.vstack ((C2, row) )
    # set up the picos problem
    prob = pic.Problem()
    K_value = K
    Q_value = Q
    # define variables
    N = Constant("N", range(len(q2))) # clients {1,2,...N} plus destination, column vector
    K = Constant("K", K_value) # number of vehicles
    Q = Constant("Q", Q_value) # vehicle capacity
    q = Constant("q", q2) # demand for each client
    c = Constant("c", C3, C3.shape) # travel cost from i to j
    x = BinaryVariable("x", C3.shape) # whether edge is in the tour
    u = RealVariable("u", len(q2)) # cumulative quantity of good delivered between origin and client i
    prob.set_objective("min", pic.sum([c[i,j]*x[i,j] for j in range(C3.shape[1]) for i in range(C3.shape[0])])) # double check this
    # add constraints
    # all u values should be >= q and <= Q
    for i in range(len(q2)):
        prob.add_constraint(u[i] >= q[i])
        prob.add_constraint(u[i] <= Q)
    # for all edges, u_i-u_j+Q*x <= Q - qj
    for i in range(C3.shape[0]) :
        for j in range(C3.shape[1]):
            prob.add_constraint(u[i] - u[j] + Q * x[i,j] <= Q - q[j])
    # make sure that only 1 vehicle leaves every client
    # make sure rows of x sum to 1
    for n in N[1:-1]:
        prob.add_constraint(pic.sum(x[n,:])==1)
    # make sure that only 1 vehicle enters every client
    # make sure cols of x sum to 1
    for n in N[1:-1]:
        prob.add_constraint(pic.sum(x[:,n])==1)
    # make sure that no more than K vehicles leave the origin
    prob.add_constraint(pic.sum(x[0,:]) <= K)
    # make sure that no more than K vehicles enter the origin
    prob.add_constraint(pic.sum(x[:,-1]) <= K)
    # make sure that the number of vehicles that leave the origin should
    # equal the number of vehicles coming back in
    prob.add_constraint(pic.sum(x[0,:]) == pic.sum(x[:,-1]))
    # solve integer program
    res = prob.solve(solver='cplex')
    objective_value = pic.sum([c[i,j]*x[i,j] for j in range(C3.shape[1]) for i in range(C3.shape[0])])
    return objective_value, x

# Local search approach (OPTIONAL)
def local_search(C,q,K,Q):
    '''
    Solves the capacitated vehicle routing problem using a local search
    approach.

    C: matrix of edge costs, that represent distances between each node
    q: list of demands associated with each client node
    K: number of vehicles
    Q: capacity of each vehicle
    returns:
        bestval: value of the minimum travel cost
        bestx: matrix representing number of routes that use each arc
    '''
    bestx = []
    bestval = 0

    # TODO (OPTIONAL): implement local search to solve vehicle routing problem

    return bestval, bestx


if __name__ == "__main__":

    # an example call to test your integer programming implementation
    C,q,K,Q = read_file_type_A('data/A-n05-k04.xml')
    travel_cost, x = cvrp_ip(C,q,K,Q)
    print("Travel cost: " + str(travel_cost))
