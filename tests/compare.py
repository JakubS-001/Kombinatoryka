#!/usr/bin/env python3
import sys, json

def read_floats(path):
    with open(path) as f:
        vals = [float(x) for x in f.read().strip().split() if x.strip()!='']
    return sorted(vals)

if __name__=='__main__':
    if len(sys.argv) < 4:
        print('Usage: compare.py target.txt computed.txt out.json [eps]', file=sys.stderr)
        sys.exit(2)
    target = read_floats(sys.argv[1])
    comp = read_floats(sys.argv[2])
    out = sys.argv[3]
    eps = float(sys.argv[4]) if len(sys.argv)>4 else 1e-6
    res = {'passed': False, 'maxdiff': None, 'n_target': len(target), 'n_computed': len(comp)}
    if len(target) != len(comp):
        res['reason'] = 'length_mismatch'
    else:
        diffs = [abs(a-b) for a,b in zip(target,comp)]
        maxd = max(diffs) if diffs else 0.0
        res['maxdiff'] = maxd
        res['passed'] = maxd <= eps
    with open(out,'w') as f:
        json.dump(res, f, indent=2)
    print(json.dumps(res))
