def ExtractContour(cont):
    p = cont.collections[0].get_paths()[0]
    v = p.vertices
    x = v[:,0]
    y = v[:,1]

    return x, y
