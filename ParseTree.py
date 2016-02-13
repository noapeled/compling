"""
Provides the ParseTree and ParseTreeNode classes for representing PCFG parse trees.

@author: Noa Peled
"""

#================================================================================
# Classes
#================================================================================

class ParseTree:
    """
    Represents a parse of a sentence according to a particular PCFG.
    
    @ivar root: The root node in the parse tree. This is usually a node for "S".
    @type root: ParseTreeNode
    @ivar probability: The probability that the PCFG assigns to this parse tree. It is the product of the 
        probabilities of all the rules used to derive this tree.
    @type probability: float 
    """
    def __init__(self, root, probability = 0):
        self.root = root
        self.probability = probability
    
    def __repr__(self):
        return "(%f): %r" % (self.probability, self.root)


class ParseTreeNode:
    """
    Represents a single node in a PCFG parse tree.
    
    @ivar key: The variable or terminal at this node.
    @type key: string
    @ivar children: The child nodes of self in the parse tree. If self is a terminal, this 
        sequence is empty. 
    @type children: A sequence (e.g. a list) of ParseTreeNode instances. 
    """
    def __init__(self, key, children = None):
        self.key = key
        
        if children == None:
            children = []
        self.children = children
    
    def __repr__(self):
        if len(self.children) > 0:
            return "[%s %s]" % (self.key, " ".join(map(repr, self.children)))
        else:
            return self.key
