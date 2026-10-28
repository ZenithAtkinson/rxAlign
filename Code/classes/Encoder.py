import jsonpickle

class Encoder:
    """
    Pickles and saves objects as json files or loads and unpickles json files into objects
    """
    def save(self, obj, file_path):
        """
        Pickles obj and saves as json at file_path.
        """
        json_string = jsonpickle.encode(obj, unpicklable=True)
        # json_data = json.dumps(json_string, indent=4)
        with open(file_path, "w") as output:
            output.write(json_string)

    def load(self, file_path):
        """
        Reads json data from file_path, unpickles it, and returns the obj 
        """
        json_data = ""
        with open(file_path, "r") as input:
            json_data = input.read()
        obj = jsonpickle.decode(json_data)
        # obj = json.loads(json_data)
        return obj
        