"""
Provides the PCFG, NearCNF and PCFGRule classes for representing and parsing probabilistic
context-free grammars.

@author: Noa Peled
"""

#================================================================================
# Imports
#================================================================================
import copy
from collections import defaultdict
from ParseTree import ParseTree, ParseTreeNode
from itertools import combinations

#================================================================================
# Constants
#================================================================================
TERMINAL_VAR_PREFIX = 'TERMINAL_'
SHORTENED_VAR_PREFIX = 'SHORTENED_'
EPSILON = ''

TERMINAL_BACK_POINTER = 'TERMINAL_BACK_POINTER'
ORDINARY_BACK_POINTER = 'ORDINARY_BACK_POINTER'
UNIT_RULE_BACK_POINTER = 'UNIT_RULE_BACK_POINTER'

#================================================================================
# Classes
#================================================================================

class PCFGBase:
    """
    An abstract base class for PCFG classes. (Do not use this class directly, use the PCFG class.) 
    
    @ivar start_variable: The variable from which all strings in the language are derived. This is typically "S".
    @type start_variable: string
    @ivar rules: The rules of the grammar.
    @type rules: A sequence of PCFGRule objects.
    """
    
    def __init__(self, start_variable="S", rules=None):
        """
        Initializes a PCFG instance with the specified set of rules and start variable.
        
        @keyword start_variable: The variable from which all strings in the language are derived. Default is "S".
        @type start_variable: string
        @keyword rules: The rules of the grammar. Default value is the empty list.
        @type rules: A sequence of PCFGRule objects.
        """
        self.start_variable = start_variable
        
        if rules is None:
            rules = []
        self.rules = rules
    
    def is_valid(self):
        """
        Checks that the grammar is legal, meaning that the probabilities of the rules for each
            variable are non-negative and sum to 1.
        
        @return: Whether the grammar is valid.
        @rtype: bool
        """
        sum_prob_per_var = {}
        for rule in self.rules:
            var, prob = rule.variable, rule.probability
            if prob < 0:
                return False
            sum_prob_per_var[var] = sum_prob_per_var.get(var, 0) + prob
        return all(sum_prob == 1.0 for sum_prob in sum_prob_per_var.values())

    @staticmethod
    def is_variable(item):
        """
        A utility function, returns True if the given item is a variable, False Otherwise.

        @param item: The item to check.
        @type item: string

        @return: Check result.
        @rtype: bool
        """
        return len(item) > 0 and item[0].capitalize() == item[0]

class PCFG(PCFGBase):
    """
    Represents a probabilistic context-free grammar.
    
    @ivar cnf: An equivalent grammar in near-CNF.
    @type cnf: NearCNF
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cnf = None
    
    def parse(self, sentence):
        """
        Parses the given sentence.
        
        @param sentence: A sequence of natural language words.
        @type sentence: A sequence (e.g. a tuple) of strings.

        @return: If sentence is in the language, returns a valid parse ParseTree of the sentence that has maximum
            probability. Otherwise, returns None.
        @rtype: ParseTree or None.
        
        @precondition: self.cnf is None or self.cnf is up-to-date. If the grammar has been changed since self.cnf 
            was computed last time, the user must set self.cnf back to None before calling parse.
        """
        if self.cnf is None:
            self.compute_cnf()
        cnf_tree = self.cnf.cky_parse(sentence)
        original_tree = self.get_original_tree(cnf_tree)
        return original_tree

# ---------------- FUNCTIONS FOR CONVERTING A GENERAL PCFG TO A NEAR-CNF EQUIVALENT -----------

    @staticmethod
    def __remove_from_rhs(rule, A, epsilon_prob, seen_epsilon_vars):
        """
        Returns all possible rules resulting from removing at least one A from the right-hand side of rule.

        @param rule: The rule where A appears in the right-hand side.
        @type rule: C{PCFGRule}

        @param A: the variable to remove from the right-hand side.
        @type A: L{string}

        @param epsilon_prob: the probability which epsilon rule A -> epsilon had.
        @type epsilon_prob: float

        @param seen_epsilon_vars: Variables from epsilon rules which we've already removed.
        @type seen_epsilon_vars: set of string

        @return: All resulting rules.
        @rtype: list of PCFGRule
        """
        derivation = rule.derivation
        if A not in derivation:
            return [rule.copy()]

        variable = rule.variable
        without_A = []
        # Define p, q as in the formula from class
        p, q = epsilon_prob, rule.probability
        indexes_with_A = [i for i in range(len(derivation)) if derivation[i] == A]

        # Do not remove all A's from right-hand side if this results in an epsilon rule which we've previously removed.
        num_iterations = len(indexes_with_A)
        if not ((variable in seen_epsilon_vars) and (all(rhs_var == A for rhs_var in derivation))):
            num_iterations += 1

        for num_indexes_to_remove in range(num_iterations):
            for indexes_to_remove in combinations(indexes_with_A, num_indexes_to_remove):
                # Define K and L as in the formula from class
                K = num_indexes_to_remove  # Number of A's removed
                L = len(indexes_with_A) - K  # Number of A's not removed
                rhs_without_A = []
                for i in range(len(derivation)):
                    if derivation[i] != A or i not in indexes_to_remove:
                        rhs_without_A.append(derivation[i])
                probability = q * (p ** K) * ((1.0 - p) ** L)
                without_A.append(PCFGRule(variable, rhs_without_A, probability,
                        {"rule": rule, "indexes_to_remove": indexes_to_remove, "removed_variable": A}))
        return without_A

    @staticmethod
    def __remove_epsilon_rule(epsilon_rule, rules, seen_epsilon_vars):
        """
        For internal use: removes a given epsilon rule from rules, per the conversion algorithm.

        @param epsilon_rule: The epsilon rule to remove.
        @type epsilon_rule: C{PCFGRule}

        @param rules: Rules to remove from.
        @type rules: list of C{PCFGRule}

        @param seen_epsilon_vars: Variables from epsilon rules which we've already removed.
        @type seen_epsilon_vars: set of string

        @return: The resulting rules
        @rtype: list of C{PCFGRule}
        """
        A = epsilon_rule.variable
        factor = 1.0 / (1.0 - epsilon_rule.probability)
        new_rules_2b = []
        for rule in rules:
            # 2.a) Don't take the given epsilon rule into the result rules.
            if rule.variable == A and not rule.derivation:
                continue
            # 2.b) Normalize probabilities of other rules.
            variable = rule.variable
            new_rule = rule.copy()
            if rule.variable == A:
                new_rule.probability *= factor
            new_rules_2b.append(new_rule)

        # 2.c) Adapt rules where A appears on the right-hand side.
        new_rules_2c = []
        for rule in new_rules_2b:
            new_rules_2c.extend(PCFG.__remove_from_rhs(rule, A, epsilon_rule.probability, seen_epsilon_vars))

        # Finally, consolidate every set of rules {X -> u [q_1], ..., X -> u [q_n]}
        # into one rule X -> u [q_1 + ... + q_n]
        repeated_rules = defaultdict(float)
        for rule in new_rules_2c:
            repeated_rules[rule.variable, tuple(rule.derivation)] += rule.probability
        new_rules = [PCFGRule(v, d, repeated_rules[v, d]) for v, d in repeated_rules]
        return new_rules

    @staticmethod
    def __find_unseen_epsilon_rule(rules, seen_epsilon_vars):
        """
        Returns a variable X, for which there exists an epsilon rule not previously removed.

        @param rules: Rules to look for X in.
        @type rules: list of PCFGRule

        @param seen_epsilon_vars: Variables for which epsilon rules have previously been removed.
        @type seen_epsilon_vars: set of string

        @return: Variable X if found, None otherwise.
        @rtype: string or NoneType
        """
        for rule in rules:
            if rule.variable not in seen_epsilon_vars and not rule.derivation:
                return rule

    def __conversion_step_2(self):
        """
        For internal use: step 2 of the conversion algorithm, removes epsilon rules and changes other rules accordingly.
        We assume there is at most one rule of the form X -> epsilon for each variable X.

        @return: the new rules,  which result from step 2 of the conversion algorithm.
        @rtype: list of C{CFGRule}
        """
        # Set of variables for which epsilon rules have been removed already
        seen_epsilon_vars = set()
        new_rules = self.rules
        while True:
            # Find a variable A which has a rule A -> epsilon and which we haven't removed already
            epsilon_rule = self.__find_unseen_epsilon_rule(new_rules, seen_epsilon_vars)
            if not epsilon_rule:
                break
            new_rules = PCFG.__remove_epsilon_rule(epsilon_rule, new_rules, seen_epsilon_vars)
            seen_epsilon_vars.add(epsilon_rule.A)
        return new_rules

    def __shorten_rule(rule):
        """
        Implements step 3 for a single rule, that is converts a rule X->u1,...,un [q] to a list
        [X-> u1X1[q] , X1 -> u2X2 [1], ..., X(n-2) -> u(n-1)un [1]].
        If the right hand side of the rule has at most two items, a list with the same rule is returned
        @param rule: the rule to convert
        @type rule: PCFGRule
        @return: the resulting rules
        @rtype: list of PCFGRule
        """
        if len(rule.derivation) <= 2:
            return [rule.copy()]

        PCFGRuleLst = []
        LHS = rule.variable
        for i in range(len(rule.derivation) - 2):
            new_var = SHORTENED_VAR_PREFIX + rule.variable + str(i + 1)
            RHS = [rule.derivation[i], new_var]
            probability = rule.probability if i == 0 else 1.0
            PCFGRuleLst.append(PCFGRule(LHS, RHS, probability, {"rule": rule}))
            LHS = new_var

        PCFGRuleLst.append(PCFGRule(LHS, rule.derivation[-2:]), 1.0, {"rule": rule})

        return PCFGRuleLst

    @staticmethod
    def __is_terminal(item):
        # Returns a boolean indicating if the item is a terminal
        return item[0] == item[0].lower()

    def __conversion_step_4(self):
        # 4a) Shorten long rules
        short_rules = []
        for rule in self.rules:
            short_rules.extend(PCFG.__shorten_rule(rule))

        # 4b) Replace terminals on right hand sides of length 2.
        new_rules = []
        terminals = set()  # Will keep track of terminals which should correspond to new variables
        terminal_to_var = lambda terminal: TERMINAL_VAR_PREFIX + item.capitalize()
        for rule in short_rules:
            rule_copy = rule.copy()
            if (len(rule_copy.derivation >= 2)):
                for i in range(2):
                    item = rule_copy.derivation[i]
                    if PCFG.__is_terminal(item):
                        terminals.add(item)
                        rule_copy.derivation[i] = terminal_to_var(item)
            new_rules.append(rule_copy)
        new_rules.extend(PCFGRule(terminal_to_var(terminal), terminal, 1.0) for terminal in terminals)
        return new_rules

    def compute_cnf(self):
        """
        Creates an equivalent near-CNF grammar and stores it in self.cnf.
        
        The equivalent grammar satisfies the following conditions:
            * It is in near-CNF.
            * It generates the same language as this PCFG (self).
            * For every string s in the language, there is a bijection between the valid parse trees of s in 
                the original grammar and the the valid parse trees of s in the new grammar such that matching  
                trees have the same probability.
                
        In addition, this method stores the information needed for converting trees from the near-CNF grammar back
        to the equivalent trees in the original grammar. This information will then be used by get_original_tree. 
        You are free to choose how this information is represented and where it is stored.
        
        @postcondition: self.cnf contains a NearCNF object with an equivalent grammar.
        """
        # Initialize an empty near-CNF, then fill it up with rules per the conversion algorithm.
        start_variable = self.start_variable + '_0'
        near = NearCNF(start_variable)

        # 1) Add rule S_0 -> S [1]
        near.rules.append(PCFGRule(start_variable, self.start_variable, 1.0))

        # 2) Remove epsilon rules. We assume there is at most one rule of the form X -> epsilon for each variable X.
        near.rules.extend(self.__conversion_step_2)

        # 3) Remove unit rules -- skipped, since unit rules are allowed in a near-CNF grammar.

        # 4) Shorten long rules and replace terminals on right-hand sides of length two.
        near.rules.extend(self.__conversion_step_4())

        return near

    # --------- FUNCTIONS FOR REVERTING BACK TO THE ORIGINAL GRAMMAR ---------

    def __revert_terminal_variables(self, root):
        """
        Reverts terminal variables in the tree rooted at node.
        @param root: The root of the tree.
        @type root: C{ParseTreeNode}
        """
        var_to_terminal = lambda var_name: var_name[len(TERMINAL_VAR_PREFIX):].lower()
        for node in root:
            if node.key.startswith(TERMINAL_VAR_PREFIX):
                node.key = var_to_terminal(node.key)
                node.children = []

    def __revert_short_rules(self, root):
        """
        Reverts a succession of short rules created during conversion back to the long rule which yielded them.
        @param root: The root of tree in which to revert the rules.
        @type root: C{ParseTreeNode}
        """
        for node in root:
            if any(child.key.startswith(SHORTENED_VAR_PREFIX) for child in node.children):
                new_children = []
                curr_node = node
                while curr_node.children[1].key.startswith(SHORTENED_VAR_PREFIX):
                    new_children.append(curr_node.children[0])
                    curr_node = curr_node.children[1]
                new_children.extend(curr_node.children)
                node.children = new_children

    def __shorten_path(self, node):
        replacement_children = []
        while True:
            for rule in self.rules:
                if rule.variable == node.key:
                    replacement_children.append(rule.children[0])
                    node = rule.children[1]
                    break
            else:
                break
        replacement_children.extend(node.children)
        return replacement_children

    def __revert_step_4(self, root):
        """
        Reverts step 4: terminal variables and short rules.

        @param root: root of the tree in which to revert.
        @type root: C{ParseTreeNode}
        """
        self.__revert_terminal_variables(root)
        self.__revert_short_rules(root)

    def __revert_step_2(self, root):
        """
        Reverts step 2: epsilon rules.

        @param root: root of the tree in which to revert.
        @type root: C{ParseTreeNode}
        """
        for node in root:
            try:
                indexes_to_revert = node.rule.original_rule.indexes_to_remove
                variable_to_revert = node.rule.original_rule.removed_variable
                for index in indexes_to_revert:
                    node.children.insert(index, ParseTreeNode(
                            variable_to_revert, [ParseTreeNode(EPSILON)]))
                node.rule.original_rule = None
            except AttributeError:
                continue

    def get_original_tree(self, tree):
        """
        Takes a valid parse ParseTree in the grammar of self.cnf and returns the equivalent parse ParseTree in 
        the the original grammar (the grammar of self).
        
        In particular, the output parse ParseTree has the same probability as the input parse ParseTree. 
        
        @param tree: The near-CNF parse ParseTree, or None.
        @type tree: ParseTree, or None.
        
        @return: If ParseTree is not None, returns the original ParseTree. Otherwise, returns None.
        @rtype: ParseTree or None.
        """
        if not tree:
            return
        tree = copy.deepcopy(tree)
        self.__revert_step_4(tree.root)
        self.__revert_step_2(tree.root)
        # Get rid of step 1, namely get rid of S_0 -> S
        new_root = tree.root.children[0]
        new_tree = ParseTree(new_root, tree.probability)
        return new_tree

class NearCNF(PCFGBase):
    """
    Represents a PCFG in near-CNF.
    """

    class __UnitRulesGraph:
        """

        """
        def __init__(self, near_cnf_grammar):
            """

            @param near_cnf_grammar:
            @return:
            """
            unit_rules = tuple(rule for rule in near_cnf_grammar.rules
                               if len(rule.derivation) == 1 and PCFG.is_variable(rule.derivation[0]))
            self.vertices = frozenset(rule.variable for rule in unit_rules)
            self.neighbors = {v: [] for v in self.vertices}
            for rule in unit_rules:
                self.neighbors[rule.variable].append(rule.derivation[0], rule.probability)

    @staticmethod
    def __dijkstra_max_prob(unit_rules_graph, source_var):
        """

        @param unit_rules_graph:
        @param source_var:
        @return:
        """
        variables = unit_rules_graph.vertices
        dist = {var: None for var in variables}
        dist[source_var] = 1.0
        tree_nodes = {var: ParseTreeNode(var) for var in variables}
        queue = set(variables)
        while queue:
            max_dist, max_dist_var = max((d, v) for v, d in dist.items() if d is not None)
            queue.remove(max_dist_var)
            for neighbor, prob in unit_rules_graph.neighbors[max_dist_var]:
                if neighbor not in queue:
                    continue
                alt_dist = dist[max_dist_var] * prob
                if dist[neighbor] is None or alt_dist > dist[neighbor]:
                    dist[neighbor] = alt_dist
                    tree_nodes[max_dist_var].children.append(tree_nodes[neighbor])
        return ParseTree(tree_nodes[source_var])

    def __compute_unit_routes(self):
        unit_rules_graph = NearCNF.__UnitRulesGraph(self)
        unit_routes = {var: None for var in unit_rules_graph.vertices}
        return unit_routes

    def __searchable_rules(self):
        searchable_rules = {rule.var: {} for rule in self.rules}
        for rule in self.rules:
            searchable_rules[rule.var][tuple(rule.derivation)] = rule.probability
        return searchable_rules

    def __best_units_derivation(searchable_rules, unit_routes, lhs_var, final_rhs):
        best_route = None
        best_route_prob = searchable_rules[lhs_var].get(final_rhs, 0.0)
        for route in unit_routes[lhs_var]:
            full_route = [lhs_var] + route
            curr_prob = 1.0
            for i in range(1, len(full_route)):
                prev_var = full_route[i - 1]
                next_var = full_route[i]
                curr_prob *= searchable_rules[prev_var][(next_var,)]
                alt_prob = curr_prob * searchable_rules[next_var].get(final_rhs, 0.0)
                if alt_prob > best_route_prob:
                    best_route = full_route[:i + 1]
                    best_route_prob = alt_prob
        return best_route, best_route_prob

    @staticmethod
    def __recursive_backtrack(back, i, j, item):
        """

        @param back:
        @param i:
        @param j:
        @param item:
        @return:
        """
        # Base case: item is a terminal.
        if not PCFG.is_variable(item):
            return ParseTreeNode(item)

        # Recursive cases: item is a variable.
        backpointer = back[i][j][item]
        children = []
        curr_children = children

        # backpointer may begin with a list of variables [A1, ..., Ar], corresponding to a chain of unit rules
        # item -> A1 -> ... -> Ar.
        if backpointer.type == UNIT_RULE_BACK_POINTER:
            chain = backpointer[0]
            # Convert the chain of unit rules to corresponding child nodes.
            for var in chain:
                next_node = ParseTreeNode(var)
                curr_children.append(next_node)
                curr_children = next_node.children
            backpointer = backpointer[1:]

        # backpointer will now have one of the following forms:
        # Case 0: (str,) = (terminal,) corresponds to directly deriving the terminal.
        if isinstance(backpointer[0], str):
            curr_children.append(NearCNF.__recursive_backtrack(back, i, j, backpointer[0]))
        # Case 1: (int, variable1, variable2) = (k, B, C) as in the CKY algorithm for CNF.
        else:
            k, B, C = backpointer
            curr_children.extend([NearCNF.__recursive_backtrack(back, i, k, B),
                                  NearCNF.__recursive_backtrack(back, k, j, C)])

        # Finally, return the node rooted at the given item.
        return ParseTreeNode(item, children)

    def __reconstruct_tree(self, probs, back, sentence_len):
        """
        Helper function for cky_parse, reconstructs the maximum probability tree from the given tables.

        @param probs: The table ("t") of probabilities computed in cky_parse.
        @type probs: list

        @param back: The table of backtrack pointers computed in cky_parse.
        @type back: list

        @param sentence_len: number of tokens ("T") in the sentence which cky_parse worked on.
        @type sentence_len: int

        @return: The reconstructed tree, if the sentence which cky_parse worked on is in the language, None otherwise.
        @rtype: C{ParseTree} or C{NoneType}
        """
        # If sentence isn't in the language of the grammar.
        if not back[0][sentence_len][self.start_variable]:
            return
        # Otherwise, reconstruct and return the computed maximum probability tree for the sentence.
        root = self.__recursive_backtrack(back, 0, sentence_len, self.start_variable)
        # The tree has probability as computed in cky_parse.
        tree_prob = probs[0][sentence_len][self.start_variable]
        return ParseTree(root, tree_prob)

    def cky_parse(self, sentence):
        """
        Parses the given text using a variant of the CKY algorithm for near-CNF grammars.

        COMPLEXITY: O(n^2*G^2 + n^3*G), where: n = #tokens in sentence, G = #rules in grammar.

        @param sentence: A sequence of natural language words.
        @type sentence: A sequence (e.g. a tuple) of strings.
        
        @return: If sentence is in the language, returns a valid parse ParseTree of the sentence that has maximum 
            probability. Otherwise, returns None.
        @rtype: ParseTree or None.
        """
        # This code is based on the variant of CKY from HW9, which can also deal with unit productions.
        # After filling a cell with variables as per the original CKY algorithm, the variant adds to the cell
        # every variable var1 such that \exists var2 in the cell so that var1 =>* var2.
        sentence = sentence.split()
        T = len(sentence)
        unit_rules, terminal_rules = self.__preprocess_unit_rules()

        # The 3D tables of dimensions (T+1)x(T+1)x|V| are each implemented as a nested list,
        # such that each cell [i][j] holds a dict which maps variables to probabilities (table t)
        # or to backtrack pointers (table back).
        t = [[{} for j in range(T + 1)] for i in range(T + 1)]
        back = [[None for j in range(T + 1)] for i in range(T + 1)]
        unit_routes = self.bfs_unit_routes()

        # Build tables.
        for j in range(1, T + 1):
            # pseudo code: "t[j - 1, j, A] = Prob(A -> Oj) for each A -> Oj in G"
            word_j = sentence[j]
            for rule in self.rules:
                if rule.derivation == [word_j]:
                    t[j - 1][j][rule.variable] = rule.probability
                    back[j - 1][j][rule.variable] = {"type": TERMINAL_BACK_POINTER, "rule": rule}

            # Improve using unit rules
            while True:
                improved = False
                for var, rules in unit_rules.items():
                    for rule in rules:
                        v = rule.variable
                        alt_prob = rule.probability * terminal_rules[v][]
                if not improved:
                    break

            for i in range(j - 2, -1, -1):
                # First, handle rules of the form A -> BC.
                for k in range(i + 1, j):
                    for rule in filter(lambda r: len(r.derivation) == 2, self.rules):
                        A = rule.variable
                        B, C = rule.derivation
                        alt_prob = rule.probability * t[i][k][B] * t[k][j][C]
                        if t[i][j][A] < alt_prob:
                            t[i][j][A] = alt_prob
                            back[i][j][A] = {"type": ORDINARY_BACK_POINTER, "k": k, "rule": rule}
                # Then, handle unit rules.
                for A, routes in unit_routes.items():
                    # Note that from here on to the end of the loop, each variable is traversed at most once,
                    # hence each iteration takes time O(G).
                    for route in routes:
                        alt_prob = 1.0
                        for var, prob in route:
                            alt_prob *= prob * t[i][j][var]
                            if t[i][j][A] < alt_prob:
                                t[i][j][A] = alt_prob
                                back[i][j][A] = {"type": UNIT_RULE_BACK_POINTER, "rule": rule}

        return self.__reconstruct_tree(t, back, T)

class PCFGRule:
    """
    Represents a single rule in a PCFG.
    
    Rules are of the form:
        variable -> derivation (probability)
    
    @ivar variable: The variable on the left side of the rule.
    @type variable: string
    
    @ivar derivation: The variable on the right side of the rule.
    @type derivation: Sequence of strings (e.g. a tuple of strings).
    
    @ivar probability: The probability that the variable will be transformed to the derivation.
    @type probability: float 
    """
    
    def __init__(self, variable, derivation, probability, original_rule=None):
        self.variable = variable
        self.derivation = derivation
        self.probability = probability
        # For rules created through conversion from other rules.
        self.original_rule = original_rule

    def copy(self):
        new_rule = copy.deepcopy(self)
        new_rule.original_rule = {"rule": self}
        return new_rule