from .expression import *

__version__ = "1.0"
__all__ = ['is_unary_operator', 'is_binary_operator',
           'is_boolean_operator', 'is_valid_variable_name',
           'BooleanParseTree', 'Stack', 'build_parse_tree']


def is_boolean_operator( token ):
    return isinstance(token, str) and (token in supported_boolean_operators())


def is_unary_operator(token):
    return isinstance(token, str) and (token in unary_boolean_operators())


def is_binary_operator(token):
    return isinstance(token, str) and (token in binary_boolean_operators())


def is_valid_variable_name(token):
    return isinstance(token, str) \
           and len(token) > 0 \
           and (token not in supported_boolean_operators()) \
           and (token.lower()[0] in "abcdefghijklmnopqrstuvwxyz")


class Stack:
    def __init__(self, max_size=None):
        self.__stack = list()
        self.__max_size = max_size

    def push(self, data):
        if (self.__max_size is None) or (len(self.__stack) < self.__max_size):
            self.__stack.append(data)

    def pop(self):
        return self.__stack.pop(-1)

    def __len__(self):
        return len(self.__stack)

    def __getitem__(self, item):
        return self.__stack[item]

    def is_empty(self):
        return len(self.__stack) == 0

    def is_full(self):
        return (self.__max_size is not None) and (len(self.__stack) < self.__max_size)


class BooleanParseTree:
    def __init__(self, parent=None):
        self.__root = None
        self.__root_type = None
        self.__parent = parent
        self.__leftChild = None
        self.__rightChild = None

    @property
    def is_top(self):
        return self.__parent is None

    @property
    def root_name(self):
        return self.__root

    def set_root(self, val):
        assert self.__root is None
        assert is_boolean_operator(val) or is_valid_variable_name(val)
        if is_boolean_operator(val):
            self.__root_type = 'Operator'
        else:
            self.__root_type = 'Variable'
        self.__root = val

    def create_left(self):
        self.__leftChild = BooleanParseTree(parent=self)

    def create_right(self):
        self.__rightChild = BooleanParseTree(parent=self)

    @property
    def has_child(self):
        return (self.__leftChild is not None) or (self.__rightChild is not None)

    @property
    def root_type(self):
        return self.__root_type

    @property
    def right_child(self):
        return self.__rightChild

    @property
    def left_child(self):
        return self.__leftChild

    @property
    def parent(self):
        return self.__parent

    @property
    def n_child(self):
        return int(self.__leftChild is not None) + int(self.__rightChild is not None)

    def is_valid(self):
        if self.__root is None:
            sts = self.n_child == 0
        elif is_boolean_operator(self.__root):
            valid_left, valid_right = False, False
            if self.__leftChild is not None:
                valid_left = self.__leftChild.is_valid()
            if self.__rightChild is not None:
                valid_right = self.__rightChild.is_valid()

            if is_unary_operator(self.__root):
                sts = valid_right
            else:
                sts = valid_left and valid_right
        else:
            sts = is_valid_variable_name(self.__root) and (self.n_child == 0)
        return sts

    @property
    def leaves(self):
        assert self.is_valid()
        if self.__root is None:
            return []
        elif (self.__leftChild is None) and (self.__rightChild is None):
            return [self]
        else:
            return self.__leftChild.leaves + self.__rightChild.leaves

    def __str__(self):
        assert self.is_valid()
        ret = ''
        if self.n_child == 0:
            ret = '%s' % self.__root
        elif is_unary_operator(self.__root):
            ret = '(%s %s)' % (self.__root, self.__rightChild)
        elif is_binary_operator(self.__root):
            ret = '(%s %s %s)' % (self.__leftChild, self.__root, self.__rightChild)
        return ret


def build_parse_tree(fpexp):
    fpexp = fpexp.strip()
    e_tree = BooleanParseTree()
    current_tree = e_tree
    token = ''
    for i in fpexp:
        if i == '(':
            current_tree.create_left()
            current_tree = current_tree.left_child
            token = ''
        elif i == ')':
            if token != '':
                current_tree.set_root(token)
            assert not current_tree.is_top
            current_tree = current_tree.parent
            token = ''
        elif i in {" ", "\t"}:
            if token != '':
                if is_boolean_operator(token):
                    current_tree = current_tree.parent
                    current_tree.set_root(token)
                    current_tree.create_right()
                    current_tree = current_tree.right_child
                elif is_valid_variable_name(token):
                    current_tree.set_root(token)
            token = ''
        else:
            token = token + i
    if token != "":
        current_tree.set_root(token)
    assert current_tree.is_top
    return e_tree

