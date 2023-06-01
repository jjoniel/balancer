"""
Chemical Equation Balancer v1.0

Created by Joniel Augustine Jerome
12/17/2022

The program parses an inputed unbalanced chemical equation lacking coefficients and 
uses a* search to assign the correct coefficients to each reactant and product.

Feats:
- highly efficient search algorithm has solved any presented equation in under 2 seconds
- most equations are solved instantaneously

Input Flexibility:
- extraneous spaces between compounds are ignored
- 8 different equation separators are supported: '===>', '--->', '==>', '-->', '=>', '->', '=', '→'
- input may be provided as a command line argument or directly

Limitations:
- the separator for the equation must be one of the 8 currently supported types
- no current support for physical state labels such as (g), (s), (l), and (aq)
- no current support for input equations with prexisting coefficients or partially filled equations

"""

from heapq import heapify, heappop, heappush
import sys
import time

alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
alpha_lower = alpha.lower()

def equation(eqn):
    # This function separates the provided equation into proudcts and reactants lists
    
    # First, it identifies the separator used in the equation
    if '===>' in eqn :
        ind = eqn.index('===>')
        reactants = eqn[:ind]
        products = eqn[ind+4:]
    elif '--->' in eqn :
        ind = eqn.index('--->')
        reactants = eqn[:ind]
        products = eqn[ind+4:]
    elif '==>' in eqn :
        ind = eqn.index('==>')
        reactants = eqn[:ind]
        products = eqn[ind+3:]
    elif '-->' in eqn :
        ind = eqn.index('-->')
        reactants = eqn[:ind]
        products = eqn[ind+3:]
    elif '=>' in eqn :
        ind = eqn.index('=>')
        reactants = eqn[:ind]
        products = eqn[ind+2:]
    elif '->' in eqn :
        ind = eqn.index('->')
        reactants = eqn[:ind]
        products = eqn[ind+2:]
    elif '=' in eqn :
        ind = eqn.index('=')
        reactants = eqn[:ind]
        products = eqn[ind+1:]
    elif "→" in eqn:
        ind = eqn.index('→')
        reactants = eqn[:ind]
        products = eqn[ind+1:]

    # Then it separates the equation pivoting on the separator
    reactant_list = reactants.strip().split('+')
    reactant_list = [r.strip() for r in reactant_list]
    product_list = products.strip().split('+')
    product_list = [p.strip() for p in product_list]

    for r in range(len(reactant_list)-1):
        print(" ", reactant_list[r], end = " + ")
    print(" ", reactant_list[-1],"=> ", end="")
    for p in range(len(product_list)-1):
        print(" ", product_list[p], end = " + ")
    print(" ", product_list[-1])
    
    # It returns the two separated lists
    return reactant_list, product_list

def combine(dict1, dict2):
    # This function combines two given dictionaries, adding up common values
    # This is used to sum up the atoms across all the reactants or all the products

    ret = {}
    for key in dict1:
        if key in dict2:
            ret[key] = dict1[key] + dict2[key]
        else:
            ret[key] = dict1[key]
    for key in dict2:
        if key not in ret:
            ret[key] = dict2[key]
    return ret

def add_multiplied(orig, addition, factor):
    # For compounds with parantheses and subscripts, this function adds the correct number of atoms
    # Ex. Mg3(PO4)2 has PO4 added with a factor of 2

    for key in addition:
        if key in orig:
            orig[key] += (addition[key] * factor)
        else:
            orig[key] = (addition[key] * factor)

def count_atoms(compound):
    # This function parses and counts the number of atoms for a given compound

    symbols = {}
    skip = 0
    for c in range(len(compound)):
        if c < skip:
            continue
        char = compound[c]
        if char == "(":
            # In the event of a parantheses or bracket, the function recursively counts up the atoms within the grouping, then adds it to the main count
            for i in range(c, len(compound)):
                if compound[i] == ")":
                    if i+1 < len(compound) and compound[i+1].isnumeric():
                        add_multiplied(symbols, count_atoms(compound[c + 1 : i]), int(compound[i + 1]))
                    else:
                        add_multiplied(symbols, count_atoms(compound[c + 1 : i]), 1)
            skip = i
        elif char == "[":
            # ''
            for i in range(c, len(compound)):
                if compound[i] == "]":
                    if i+1 < len(compound) and compound[i+1].isnumeric():
                        add_multiplied(symbols, count_atoms(compound[c + 1 : i]), int(compound[i + 1]))
                    else:
                        add_multiplied(symbols, count_atoms(compound[c + 1 : i]), 1)
            skip = i
        elif char in alpha:
            symb = ""
            if c + 1 < len(compound) and compound[c + 1] in alpha_lower:
                symb += compound[c : c + 2]
                c += 1 
            else:
                symb += compound[c]
            if symb in symbols:
                d = 1
                if c + d < len(compound) and compound[c + d].isnumeric():
                    while c + d + 1 < len(compound) and compound[c + d + 1].isnumeric():
                        d += 1
                    symbols[symb] += int(compound[c + 1: c + d + 1])
                else:
                    symbols[symb] += 1
            else:
                d = 1
                if c + d < len(compound) and compound[c + d].isnumeric():
                    while c + d + 1 < len(compound) and compound[c + d + 1].isnumeric():
                        d += 1
                    symbols[symb] = int(compound[c + 1: c + d + 1])
                else:
                    symbols[symb] = 1

    # It returns a dictionary of each atom within the compound and its corresponding count
    return symbols

def compound_atoms(compound):
    # Given a compound, this function identifies the atoms it is composed of
    symbols = set()
    for c in range(len(compound)):
        char = compound[c]
        if char in alpha:
            symb = ""
            if c + 1 < len(compound) and compound[c + 1] in alpha_lower:
                symb += compound[c : c + 2]
                c += 1 
            else:
                symb += compound[c]
            symbols.add(symb)
    
    # It returns a set of the unique atomic symbols
    return symbols

def score(r_list, p_list):
    # This function scores a given set of possible coefficients
    # This serves as the heuristic for the a* search

    score = 0
    for el in r_list:
        score += abs(r_list[el] - p_list[el])
    if score % 2 == 1:
        return score * 1.2
    else:
        return score

def possible_changes(r_list, p_list, r, p, r_counts, p_counts, r_atoms, p_atoms):
    # This function generates potential new coefficients for the equation based on the currently assigned coefficients
    # This serves as a child-generating function for the a* search

    possiblesR = []
    for x in range(len(r)):
        for i in r_atoms[r_list[x]]:
            if r_counts[i] > p_counts[i] and r[x] - 2 > 0:
                possiblesR.append(r[0:x] + [r[x]-2] + r[x+1:])
            elif r_counts[i] < p_counts[i]:
                possiblesR.append(r[0:x] + [r[x]+2] + r[x+1:])
            else:
                possiblesR.append(r[0:x] + [r[x]+1] + r[x+1:])
                if r[x] - 1 > 0:
                    possiblesR.append(r[0:x] + [r[x]-1] + r[x+1:])

    possiblesP = []
    for y in range(len(p)):
        for i in p_atoms[p_list[y]]:
            if r_counts[i] < p_counts[i] and p[y] - 2 > 0:
                possiblesP.append(p[0:y] + [p[y]-2] + p[y+1:])
            elif r_counts[i] > p_counts[i]:
                possiblesP.append(p[0:y] + [p[y]+2] + p[y+1:])
            else:
                possiblesP.append(p[0:y] + [p[y]+1] + p[y+1:])
                if p[y] - 1 > 0:
                    possiblesP.append(p[0:y] + [p[y]-1] + p[y+1:])

    # It returns a list of possible changes to the product and reactant coefficients
    return possiblesR, possiblesP

def count_compounds(reactant_list, product_list, reactant_coefficients, product_coefficients):
    # This function counts the atoms in the products and reactants after adding coefficients to check for balance

    reactant_counts = {}
    product_counts = {}
    x = 0
    for r in reactant_list:
        add_multiplied(reactant_counts, count_atoms(r), reactant_coefficients[x])
        x += 1
    y = 0
    for p in product_list:
        add_multiplied(product_counts, count_atoms(p), product_coefficients[y])
        y += 1
    
    # It returns a dictionary for both the reactants and the products mapping each atom to its count
    return reactant_counts, product_counts

def simplify(coeffs):
    # This function simplifies the coefficients in the final equation if they have a common divisor

    minimum = min(coeffs[0] + coeffs[1])
    for x in range(minimum, 1, -1):
        divisible = True
        for i in coeffs[0]:
            if i % x != 0:
                divisible = False
                break
        if not divisible:
            continue
        for i in coeffs[1]:
            if i % x != 0:
                divisible = False
                break
        if(divisible):
            return [z//x for z in coeffs[0]], [y//x for y in coeffs[1]]
        
    return coeffs

def main():
    # This function controls the overall flow of the program

    if len(sys.argv) > 1:
        to_solve = sys.argv[1]

    else:
        print("Enter equation to be solved [ Ex. Na3PO4+MgCl2=NaCl+Mg3(PO4)2 ]: ", end="")
        to_solve = input()

        # Sample Test Cases
        # to_solve = "C3H8+O2=H2O+CO2"
        # to_solve = "Al+HCl=AlCl3+H2"
        # to_solve = "Na3PO4+MgCl2=NaCl+Mg3(PO4)2"

    r_list, p_list = equation(to_solve)

    r_coefficients = [1 for r in r_list]
    p_coefficients = [1 for p in p_list]
    reactant_counts, product_counts = count_compounds(r_list, p_list, r_coefficients, p_coefficients)

    reactant_atoms = {r : compound_atoms(r) for r in r_list}
    product_atoms = {p : compound_atoms(p) for p in p_list}

    if product_counts.keys() != reactant_counts.keys():
        print("Cannot be balanced!")
        exit(0)
    start = time.perf_counter()
    r_coefficients, p_coefficients = balance(r_list, p_list, r_coefficients, p_coefficients, reactant_atoms, product_atoms, start)
    end = time.perf_counter()
    
    # After balancing, the balanced equation is printed in an easy to read format.
    for r in range(len(r_list)-1):
        print("\033[0;34m", end = "")
        print(r_coefficients[r], end="\033[0m ")
        print(r_list[r], end = " + ")
    print("\033[0;34m", end = "")
    print(r_coefficients[-1], end="\033[0m ")
    print(r_list[-1],"=> ", end="")
    for p in range(len(p_list)-1):
        print("\033[0;34m", end = "")
        print(p_coefficients[p], end="\033[0m ")
        print(p_list[p], end = " + ")
    print("\033[0;34m", end = "")
    print(p_coefficients[-1], end="\033[0m ")
    print(p_list[-1])
    print("Time Elapsed:", round(end-start,1), " seconds")

def balance(reactant_list, product_list, reactant_coefficients, product_coefficients, reactant_atoms, product_atoms, start):
    # This function performs the a* star search to balance the equation

    fringe = []
    heapify(fringe)
    visited = []
    p, q = count_compounds(reactant_list, product_list, reactant_coefficients, product_coefficients)
    new_node = (score(p, q), (reactant_coefficients, product_coefficients))
    heappush(fringe, new_node)
    visited.append((reactant_coefficients, product_coefficients))
    counter = 0
    curr = 0
    pre_time = time.perf_counter()
    while fringe:
        if(time.perf_counter() - pre_time >= 10):
            print(round(time.perf_counter()-start), "seconds have passed")
            pre_time = time.perf_counter()
        counter += 1
        v = heappop(fringe)
        v = v[1]
        if(counter % 100 == 0):
            for i in range(len(v[0])):
                v[0][i] = int(v[0][i] * 1.8)
            for i in range(len(v[1])):
                v[1][i] = int(v[1][i] * 1.8)
            fringe = []
            visited = [v]
        a, b = count_compounds(reactant_list, product_list, v[0], v[1])
        curr = score(a, b)
        if(curr == 0):
            v = simplify(v)
            return v
        r_poss, p_poss = possible_changes(reactant_list, product_list, v[0], v[1], a, b, reactant_atoms, product_atoms)
        for child in r_poss:
            for child2 in p_poss:
                if (child, child2) not in visited:
                    c, d = count_compounds(reactant_list, product_list, child, child2)
                    new_node = (score(c, d), (child, child2))
                    heappush(fringe, new_node)
                    visited.append((child, child2))

if __name__ == "__main__":
    main()