class DictWrap(dict):
    def __init__(self, d):
        dict.__init__(self, d)
        for k,v in d.iteritems():
            setattr(self, k, v)
