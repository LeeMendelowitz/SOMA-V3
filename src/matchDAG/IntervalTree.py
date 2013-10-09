"""
Interval Tree downloaded from:
http://forrst.com/posts/Interval_Tree_implementation_in_python-e0K
"""
from operator import attrgetter

get_begin = attrgetter('begin')
get_end = attrgetter('end')

class IntervalTree:
    def __init__(self, intervals):
        self.top_node = self.divide_intervals(intervals)

    def divide_intervals(self, intervals):

        if not intervals:
            return None

        x_center = self.center(intervals)

        s_center = []
        s_left = []
        s_right = []

        for k in intervals:
            if k.end < x_center:
                s_left.append(k)
            elif k.begin > x_center:
                s_right.append(k)
            else:
                s_center.append(k)

        return Node(x_center, s_center, self.divide_intervals(s_left), self.divide_intervals(s_right))
        

    def center(self, intervals):
        fs = sort_by_begin(intervals)
        length = len(fs)
        return fs[int(length/2)].begin

    def search(self, begin, end=None):
        if end:
            result = []

            #for j in xrange(begin, end+1):
            for j in (begin, end):
                for k in self.search(j):
                    result.append(k)
                result = list(set(result))
            return sort_by_begin(result)
        else:
            return self._search(self.top_node, begin, [])

    def _search(self, node, point, result):
        for k in node.s_center:
            if k.begin <= point <= k.end:
                result.append(k)
        if point < node.x_center and node.left_node:
            for k in self._search(node.left_node, point, []):
                result.append(k)
        if point > node.x_center and node.right_node:
            for k in self._search(node.right_node, point, []):
                result.append(k)
        return list(set(result))

class Interval:
    def __init__(self, begin, end):
        if begin > end:
            raise RuntimeError('Begin must be less than or equal to end')
        self.begin = begin
        self.end = end
        
    def __str__(self):
        return '(%f, %f)'%(self.begin, self.end)

    def __repr__(self):
        return "Interval %s"%str(self)

class Node:
    def __init__(self, x_center, s_center, left_node, right_node):
        self.x_center = x_center
        self.s_center = sort_by_begin(s_center)
        self.left_node = left_node
        self.right_node = right_node

def sort_by_begin(intervals):
    return sorted(intervals, key=get_begin)
