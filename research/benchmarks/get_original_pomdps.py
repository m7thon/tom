import urllib.request
from io import BytesIO
import gzip

# Download POMDPs from "Tony's POMDP Page"
base_url = 'http://cs.brown.edu/research/ai/pomdp/examples/'
pomdps = {'tiger': 'tiger.aaai.POMDP', 'paint': 'paint.95.POMDP', 'bridge': 'bridge-repair.POMDP.gz', 'network': 'network.POMDP',
          'shuttle': 'shuttle.95.POMDP', 'maze4x3': '4x3.95.POMDP', 'cheese': 'cheese.95.POMDP'}
for pomdp, fname in pomdps.items():
    with urllib.request.urlopen(base_url + fname) as response, open(pomdp + '.pomdp', 'wb') as f:
        if fname[-3:] == '.gz':
            f.write(gzip.decompress(response.read()))
        else:
            f.write(response.read())
