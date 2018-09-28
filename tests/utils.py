import json

def load_results(jsonfile):
    with open(jsonfile) as j:
        return json.loads(j.read())

