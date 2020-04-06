import picos as pic
from picos import RealVariable
from copy import deepcopy
from heapq import *
import heapq as hq
import numpy as np
import itertools
import math
from collections import deque
counter = itertools.count()

class BBTreeNode():
    def __init__(self, vars = [], constraints = [], objective='', prob=None):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.children = []
        self.prob = prob

    def __deepcopy__(self, memo):
        '''
        Deepcopies the picos problem
        This overrides the system's deepcopy method bc it doesn't work on classes by itself
        '''
        newprob = pic.Problem.clone(self.prob)
        return BBTreeNode(self.vars, newprob.constraints, self.objective, newprob)

    def buildProblem(self):
        '''
        Bulids the initial Picos problem
        '''
        prob=pic.Problem()

        prob.add_list_of_constraints(self.constraints)

        prob.set_objective('max', self.objective)
        self.prob = prob
        return self.prob

    def is_integral(self):
        '''
        Checks if all variables (excluding the one we're maxing) are integers
        '''
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
                return False
        return True

    def branch_floor(self, branch_var):
        '''
        Makes a child where xi <= floor(xi)
        '''
        n1 = deepcopy(self)
        n1.prob.add_constraint( branch_var <= math.floor(branch_var.value) ) # add in the new binary constraint

        return n1

    def branch_ceil(self, branch_var):
        '''
        Makes a child where xi >= ceiling(xi)
        '''
        n2 = deepcopy(self)
        n2.prob.add_constraint( branch_var >= math.ceil(branch_var.value) ) # add in the new binary constraint
        return n2


    def bbsolve(self):
        '''
        Use the branch and bound method to solve an integer program
        This function should return:
            return bestres, bestnode_vars

        where bestres = value of the maximized objective function
              bestnode_vars = the list of variables that create bestres
        '''
        # Step 1: Initialization
        # these lines build up the initial problem and adds it to the queue
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        queue = deque()
        queue.append(root)
        bestres = -1e20 # a small arbitrary initial best objective value
        bestnode_vars = root.vars # initialize bestnode_vars to the root vars

        # print('Finished initialization')

        while len(queue) > 0: # while queue isn't empty means there are still LP probs to solve
            temp_obj = queue.popleft()
            if temp_obj.is_integral():
                # print('Found potential solution')
                if temp_obj.vars[-1] > bestres:
                    # prune by optimality (found potential solution)
                    bestres = temp_obj.vars[-1]
                    bestnode_vars = temp_obj.vars
                continue
            # BRANCHING
            # choose a variable, x_j, that isn't a integer and create 2 new LP problems
            # adding constraints x_j <= floor(x_j) and x_j >= floor(x_j) + 1
            for v in temp_obj.vars[:-1]:
                if v.value is not None and abs(round(v.value) - float(v.value)) > 1e-4:
                    # print('Found non integer value:', str(v.value))
                    x_j = v
                    LP_floor_obj = temp_obj.branch_floor(x_j)
                    LP_ceiling_obj = temp_obj.branch_ceil(x_j)
                    # BOUNDING
                    # solve LP_floor
                    try:
                        floor_res = LP_floor_obj.prob.solve(solver='cvxopt')
                        if LP_floor_obj.vars[-1] > bestres: # prune by bounds
                            # print('Adding LP_floor to queue')
                            queue.append(LP_floor_obj)
                    except:
                        # prune by infeasibility
                        pass
                    # solve LP_ceiling
                    try:
                        ceiling_res = LP_ceiling_obj.prob.solve(solver='cvxopt')
                        if LP_ceiling_obj.vars[-1] > bestres: # prune by bounds
                            # print('Adding LP_ceiling to queue')
                            queue.append(LP_ceiling_obj)
                    except:
                        # prune by infeasibility
                        pass
                    break # move on to next elm in queue

        print('Final objective value:', bestres)
        return bestres, bestnode_vars
