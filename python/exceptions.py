
class BinaryParameterMissing(Exception):
    def __init__(self, parameter=None, message=None):
        if message is None:
            if parameter is not None:
                self.message = "Missing binary parameter:", parameter
            else:
                self.message = "Missing sufficient binary parameters"
        else:
            self.message = message

        super().__init__(self.message)



