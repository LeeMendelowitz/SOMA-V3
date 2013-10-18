from datetime import datetime
import sys

stderr = sys.stderr.write
flush = sys.stderr.flush

def w(msg):
    stderr(msg)
    flush()


class Clock(object):

    def __init__(self):
        self.elapsed = None
        self.start()

    def start(self):
        self.start = datetime.now()

    begin = start

    def tick(self):
        d = (datetime.now() - self.start)
        w('Elapsed: %.2f\n'%d.total_seconds())

    def end(self):
        self.end = datetime.now()
        d = (datetime.now() - self.start)
        self.elapsed = d.total_seconds()
        w('Elapsed: %.2f\n'%self.elapsed)

    stop = end
