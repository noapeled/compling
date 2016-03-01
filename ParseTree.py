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
        """
        Constructs a derivation tree.
        @param root: The root node of the tree.
        @type root: C{ParseTreeNode}
        @param probability: The probability of the derivation tree.
        @type probability: L{float}
        """
        self.root = root
        self.probability = probability
    
    def __repr__(self):
        """
        Converts a tree to string representation.
        """
        return "(%s): %r" % (self.probability, self.root)

class ParseTreeNode:
    """
    Represents a single node in a PCFG parse tree.

    @ivar key: The variable or terminal at this node.
    @type key: string

    @ivar rule: The rule in which this node is the left-hand side, if this node is internal.
    @type rule: C{PCFGRule}

    @ivar children: The child nodes of self in the parse tree. If self is a terminal, this
        sequence is empty.
    @type children: A sequence (e.g. a list) of ParseTreeNode instances.
    """
    def __init__(self, key, children=None, rule=None):
        """
        Constructs a new node. See class documentation for parameters.
        """
        self.key = key
        self.rule = rule
        if children == None:
            children = []
        self.children = children

    def __repr__(self):
        """
        Converts a node to string representation.
        """
        if len(self.children) > 0:
            return "[%s %s]" % (self.key, " ".join(map(repr, self.children)))
        else:
            return self.key

    def preorder(self):
        """
        Generator function, yields the tree nodes in pre-order traversal order, starting from the tree root.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(node.children[::-1])
