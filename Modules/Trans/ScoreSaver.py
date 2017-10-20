class ScoreSaver():
    def __init__(self, obj, method):
        self.obj = obj
        self.method = method
    
    def __call__(self, restraintId, score):
        getattr(self.obj, self.method)(restraintId, score)
