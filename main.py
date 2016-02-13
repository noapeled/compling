"""
Advance Computational Linguistics - Final Project
Tel Aviv Univerity, course no. 0627.4090, Fall 2015.

An implementation of PCFG parsing.

@author: Noa Peled
"""

#================================================================================
# Imports
#================================================================================

from PCFG import PCFG


#================================================================================
# Functions
#================================================================================

def parse(text, grammar):
    """
    Uses a PCFG grammar to parse a text.
    
    @param text: A text in natural language. You may assume that it contains only words in
        lower case separated by spaces, with periods between sentences.
    @type text: string
    
    @param garmmar: A PCFG grammar
    @type grammar: PCFG
    
    @return: Parse trees, one for each sentence in the text. The i-th tree is a parse with 
        maximum probability for the i-th sentence in the text according to the provided grammar.
    @rtype: An iterable (e.g. a list) of ParseTree instances.
    
    @todo: Implement
    """
    raise NotImplementedError()


#================================================================================
# Main
#================================================================================

if __name__ == "__main__":
    
    #
    # Todo
    # ====
    # Put a short example here. The example should initialize a PCFG instance and use it to parse
    # a short text with 5 sentences using the above parse() function. Each parse tree in the result 
    # should be printed to stdout (printed to the screen).
    #
    
    raise NotImplementedError()

