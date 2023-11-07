import numpy as np

def mat_degree(mat):
    mat = np.array(mat)
    nz = mat != 0
    return nz.sum(axis=1)

def mat_strength(mat):
    mat = np.array(mat)
    return mat.sum(axis=1)

def mat_strata(mat,n_strata=5):
    all_strata = []
    for i in range(n_strata):
        all_strata.append(np.diag(mat,i))
    return np.concatenate(all_strata)

def mat_strength_region(mat,bin,intra=True):
    result = []
    mat = np.array(mat)
    for i,j in bin[:-1],bin[1:]:
        if intra:
            average = mat[i:j,i:j].mean()
        else:
            average =  mat[i:j].mean()
        result.append(average)
    return result

def min_max_norm(mat):
    mat = (mat - mat.min())/(mat.max()-mat.min())
    return mat

def log_norm(mat):
    mat = np.log(mat + 1)
    return mat    

def mat_features(mat, operations, **kwargs):
    for op in operations:
        if op == 'degree':
            res = mat_degree(mat)
        elif op == 'strength':
            res = mat_strength(mat)
    return res

