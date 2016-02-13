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
from ParseTree import ParseTree
from itertools import combinations

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
        
        if rules == None:
            rules = []
        self.rules = rules
    
    def is_valid(self):
        """
        Checks that the grammar is legal, meaning that the probabilities of the rules for each
            variable are non-negative and sum to 1.
        
        @return: Whether the grammar is valid.
        @rtype: bool
        
        @todo: Implement.
        """
        raise NotImplementedError()


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
        if self.cnf == None:
            self.compute_cnf()
        
        cnf_tree = self.cnf.cky_parse(sentence)
        
        return self.get_original_tree(cnf_tree)

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
            return [copy.deepcopy(rule)]

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
                without_A.append(PCFGRule(variable, rhs_without_A, probability))
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
            new_rule = copy.deepcopy(rule)
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
            return [copy.deepcopy(rule)]

        PCFGRuleLst = []
        LHS = rule.variable
        for i in range(len(rule.derivation)-2):
            new_var = LHS + str(i + 1)
            RHS = [rule.derivation[i], new_var]
            probability = rule.probability if i == 0 else 1.0
            PCFGRuleLst.append(PCFGRule(LHS, RHS, probability))
            LHS = new_var

        PCFGRuleLst.append(PCFGRule(LHS, rule.derivation[-2:]), 1.0)

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
        terminal_to_var = lambda terminal: 'TERMINAL_' + item.capitalize()
        for rule in short_rules:
            rule_copy = copy.deepcopy(rule)
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
        # TODO: also keep information for mapping derivation trees.

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

    def get_original_tree(self, tree):
        """
        Takes a valid parse ParseTree in the grammar of self.cnf and returns the equivalent parse ParseTree in 
        the the original grammar (the grammar of self).
        
        In particular, the output parse ParseTree has the same probability as the input parse ParseTree. 
        
        @param ParseTree: The near-CNF parse ParseTree, or None.
        @type ParseTree: ParseTree, or None.
        
        @return: If ParseTree is not None, returns the original ParseTree. Otherwise, returns None.
        @rtype: ParseTree or None.
        
        @todo: Implement.
        """
        raise NotImplementedError()


class NearCNF(PCFGBase):
    """
    Represents a PCFG in near-CNF.
    """
    
    def cky_parse(self, sentence):
        """
        Parses the given text using a variant of the CKY algorithm for near-CNF grammars.
        
        @param sentence: A sequence of natural language words.
        @type sentence: A sequence (e.g. a tuple) of strings.
        
        @return: If sentence is in the language, returns a valid parse ParseTree of the sentence that has maximum 
            probability. Otherwise, returns None.
        @rtype: ParseTree or None.
        
        @todo: Implement.
        """
        raise NotImplementedError()


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
    
    def __init__(self, variable, derivation, probability):
        self.variable = variable
        self.derivation = derivation
        self.probability = probability