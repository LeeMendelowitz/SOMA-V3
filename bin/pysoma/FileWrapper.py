class FileWrapper(object):
    """
    Wrapper around a file object.
    If a str is passed as f, the wrapper will open and close the file.
    If a file handle is passed as f, the wrapper will not open/close the file.
    """

    def __init__(self, f, mode):
        if isinstance(f, str):
            self.file = open(f, mode)
            self.close_file = True
        else:
            self.file = f
            self.close_file = False

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        if (not self.close_file):
            return  # do nothing
        # clean up
        exit = getattr(self.file, '__exit__', None)
        if exit is not None:
            return exit(*args, **kwargs)
        else:
            exit = getattr(self.file, 'close', None)
            if exit is not None:
                exit()

    def __getattr__(self, attr):
        if (attr == 'close') and not self.close_file:
            return lambda: None
        return getattr(self.file, attr)

    def __iter__(self):
        return iter(self.file)
