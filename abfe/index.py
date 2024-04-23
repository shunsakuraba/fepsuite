

def load_index(ndx):
    group = None
    ret = {}
    with open(ndx) as fh:
        for l in fh:
            ls = l.split()
            if l.startswith('['):
                group = ls[1]
                if group in ret:
                    raise RuntimeError("Multiple groups in index file")
                else:
                    ret[group] = []
                continue
            else:
                for x in ls:
                    ret[group].append(int(x) - 1)
    return ret

