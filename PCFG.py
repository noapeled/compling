"""
Provides the PCFG, NearCNF and PCFGRule classes for representing and parsing probabilistic
context-free grammars.

@author: Noa Peled
"""

#================================================================================
# Imports
#================================================================================

from ParseTree import ParseTree


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
                    
        @todo: Implement.
        """
        raise NotImplementedError()
    
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

