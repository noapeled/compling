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
    """
    sentences = [s.rstrip() for s in text.split('.')]
    derivation_trees = [grammar.parse(sentence) for sentence in sentences]
    return derivation_trees

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
        PCFGRule('S', ['NP', 'VP'], 0.9),
        PCFGRule('S', ['VP', 'NP'], 0.1),

        PCFGRule('NP', ['DET', 'N'], 0.5),
        PCFGRule('NP', ['N'], 0.5),

        PCFGRule('VP', ['V', 'NP'], 0.3),
        PCFGRule('VP', ['V'], 0.1),
        PCFGRule('VP', ['V', 'PP'], 0.2),
        PCFGRule('VP', ['V', 'ADV'], 0.4),

        PCFGRule('PP', ['P', 'NP'], 1.0),

        PCFGRule('DET', ['the'], 0.5),
        PCFGRule('DET', ['a'], 0.25),
        PCFGRule('DET', ['an'], 0.25),

        PCFGRule('ADJ', ['few'], 0.33),
        PCFGRule('ADJ', ['old'], 0.33),
        PCFGRule('ADJ', ['prime'], 0.34),

        PCFGRule('N', ['ADJ', 'N'], 0.25),
        PCFGRule('N', ['man'], 0.125),
        PCFGRule('N', ['old'], 0.125),
        PCFGRule('N', ['boat'], 0.0625),
        PCFGRule('N', ['prime'], 0.125),
        PCFGRule('N', ['captains'], 0.125),
        PCFGRule('N', ['seven'], 0.125),
        PCFGRule('N', ['number'], 0.0625),

        PCFGRule('ADV', ['few'], 1.0),

        PCFGRule('P', ['in'], 1.0),

        PCFGRule('V', ['man'], 0.1),
        PCFGRule('V', ['is'], 0.3),
        PCFGRule('V', ['are'], 0.3),
        PCFGRule('V', ['number'], 0.1),
        PCFGRule('V', ['become'], 0.2),
    ])
    assert example_grammar.is_valid()

    trees = parse(example_text, example_grammar)
    for tree in trees:
        print(tree)
