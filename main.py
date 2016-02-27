"""
Advance Computational Linguistics - Final Project
Tel Aviv Univerity, course no. 0627.4090, Fall 2015.

An implementation of PCFG parsing.

@author: Noa Peled
"""

#================================================================================
# Imports
#================================================================================

from PCFG import PCFG, PCFGRule


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
    
    # Following is a short example text with 5 sentences and a corresponding example grammar.
    # The text is parsed per the grammar, and the resulting trees are printed.
    # Note that first and third sentences are ambiguous garden path sentences,
    # each having more than one possible derivation.
    example_text = 'the old man the boat. ' + \
                   'an old man is in the boat. ' + \
                   'the prime number few. ' + \
                   'seven is a prime number. ' + \
                   'few captains become a prime man.'

    example_grammar = PCFG(start_variable="S", rules=[
        PCFGRule('S', ['NP', 'VP']),

        PCFGRule('NP', ['DET', 'N']),
        PCFGRule('NP', ['N']),

        PCFGRule('VP', ['V', 'NP']),
        PCFGRule('VP', ['V']),
        PCFGRule('VP', ['V', 'PP']),
        PCFGRule('VP', ['V', 'ADV']),

        PCFGRule('PP', ['P', 'NP']),

        PCFGRule('DET', ['the']),
        PCFGRule('DET', ['a']),
        PCFGRule('DET', ['an']),

        PCFGRule('ADJ', ['few']),
        PCFGRule('ADJ', ['old']),
        PCFGRule('ADJ', ['prime']),

        PCFGRule('N', ['ADJ', 'N']),
        PCFGRule('N', ['man']),
        PCFGRule('N', ['old']),
        PCFGRule('N', ['boat']),
        PCFGRule('N', ['prime']),
        PCFGRule('N', ['captains']),
        PCFGRule('N', ['seven']),
        PCFGRule('N', ['number']),

        PCFGRule('ADV', ['few']),

        PCFGRule('P', ['in']),

        PCFGRule('V', ['man']),
        PCFGRule('V', ['is']),
        PCFGRule('V', ['are']),
        PCFGRule('V', ['number']),
        PCFGRule('V', ['become']),
    ])

    trees = parse(example_text, example_grammar)
    for tree in trees:
        print(tree)
