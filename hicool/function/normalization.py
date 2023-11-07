import cooler
import numpy as np
from scipy.signal import convolve2d


def fast_oe(path, chrom='chr1', method="liberman"):
    from scipy.sparse import coo_matrix
    clr = cooler.Cooler(path)
    pixels = clr.matrix(as_pixels=True).fetch(chrom)
    matrix_size = max(pixels['bin1_id'].max(), pixels['bin2_id'].max()) + 1
    # Create a sparse matrix
    matrix = coo_matrix((pixels["count"], (pixels["bin1_id"], pixels["bin2_id"])), shape=(matrix_size, matrix_size))
    # Calculate chromosome length
    chr_len = matrix_size
    # Calculate expected values for each distance
    diagonals = np.array([matrix.diagonal(i) for i in range(matrix_size)])
    if method == 'rao':
        # 对角线上的交互总和/bin数量
        expected = np.array([np.sum(diag) / len(diag) for diag in diagonals])
    elif method == 'liberman':
        # 对角线上的交互总和/染色体长度减去该距离
        expected = np.array([np.sum(diag) / (chr_len - i) for i, diag in enumerate(diagonals)])
    elif method == 'nonzero':
        # 对角线上的交互总和/该距离上的非零元素数量
        expected = np.array([np.sum(diag) / np.count_nonzero(diag) for diag in diagonals])
    elif method == 'withgap':
        pixels["distance"] = pixels["bin2_id"] - pixels["bin1_id"]
        num_vector = pixels["bin2_id"].value_counts() + pixels["bin1_id"].value_counts()
        mask = num_vector[num_vector > chr_len * 0.1]
        pixels = pixels.loc[pixels["bin1_id"].isin(mask.index) & pixels["bin2_id"].isin(mask.index)]
        d_sum = pixels.groupby("distance")["count"].sum()
        expected = d_sum / (chr_len - d_sum.index + 1)
        d_weight = d_sum / d_sum.sum()
        pixels["OE"] = (pixels["count"] * d_weight.get(pixels["distance"], 0) / expected[pixels["distance"]]).astype('float32')
        return pixels
    else:
        raise ValueError("Method must be 'rao', 'liberman', 'nonzero' or 'withgap' ")
    # Calculate OE values
    pixels["OE"] = pixels["count"] / expected[pixels["bin2_id"] - pixels["bin1_id"]]
    return pixels


def fast_norm(path, chrom='chr1', method="vc"):
    clr = cooler.Cooler(path)
    pixels = clr.matrix(as_pixels=True).fetch(chrom)
    if method == 'vc':
        count_sum = pixels.groupby(['bin1_id','bin2_id'])['count'].transform('sum')
        count_sum = count_sum.replace(0, 1)
        pixels['VC'] = pixels['count'] / (count_sum['bin1_id'] * count_sum['bin2_id'])
    elif method == 'vc_sqrt':
        count_sum = pixels.groupby(['bin1_id','bin2_id'])['count'].transform('sum')
        count_sum = count_sum.replace(0, 1)
        count_sum = np.sqrt(count_sum)
        pixels['VC_SQRT'] = pixels['count'] / (count_sum['bin1_id'] * count_sum['bin2_id'])
    return pixels



# From scHiCTools
def matrix_operation(mat, operations, **kwargs):
    """
    " Liu, J., Lin, D., Yardlmcl, G. G., & Noble, W. S. (2018). 
    Unsupervised embedding of single-cell Hi-C data. Bioinformatics "
    
    Parameters
    ----------
    mat : numpy.ndarray
        Matrix to apply smoothing operators.
        
    operations : list[str]
        A list of smoothing operators to apply.
        Now support perators: 'oe_norm', 'vc_norm', 'vc_sqrt_norm',
        'kr_norm', 'convolution', 'random_walk', 'network_enhancing',
        'logarithm', 'power'.
        
    **kwargs :
        Arguments for specific smoothing operators.
        
        'kr_norm':
            'maximum_error_rate': error rate to stop iteration, default=1e-4.
        'convolution':
            'kernel_shape': shape of kernel ('kernel_shape' ,'kernel_shape'), default=3.
        'random_walk':
            'random_walk_ratio': propotion of random walk smooth matrix, default=1.0;
            't': number of random walk iteration times, default=1.
        'network_enhancing':
            'kNN': number of nearest neighbors, default=20;
            'iterations': number of iterations, default=1.
        'logarithm':
            'epsilon': numerator of error term, default=1;
            'log_base': denominator of error's exponentiation, default=e.
        'power':
            'pow': exponent of exponentiation, default=0.5.
        

    Returns
    -------
    mat : numpy.ndarray
        Matrix after appling smoothing operator.

    """
    
    for op in operations:
        op = op.lower()
        
        if op == 'oe_norm':
            new_mat = np.zeros(mat.shape)
            averages = np.array([np.mean(mat[i:, :len(mat) - i]) for i in range(len(mat))])
            averages = np.where(averages == 0, 1, averages)
            for i in range(len(mat)):
                for j in range(len(mat)):
                    d = abs(i - j)
                    new_mat[i, j] = mat[i, j] / averages[d]
            mat = new_mat
            
        elif op == 'vc_norm':
            sm = np.sum(mat, axis=0)
            sm = np.where(sm == 0, 1, sm)
            sm_v = np.tile(sm, (len(sm), 1))
            sm_c = sm_v.T
            mat = mat / sm_c / sm_v
        
        elif op == 'vc_sqrt_norm':
            sm = np.sum(mat, axis=0)
            sm = np.where(sm == 0, 1, sm)
            sm = np.sqrt(sm)
            sm_v = np.tile(sm, (len(sm), 1))
            sm_c = sm_v.T
            mat = mat / sm_c / sm_v
            
        elif op == 'kr_norm':
            mat = KR_norm(mat, kwargs.pop('maximum_error_rate', 1e-4))
            
        elif op == 'convolution':
            mat = convolution(mat, kwargs.pop('kernel_shape', 3))
            
        elif op == 'random_walk':
            mat = random_walk(mat, kwargs.pop('random_walk_ratio', 1.0),kwargs.pop('t', 1))
            
        elif op == 'network_enhancing':
            mat = network_enhancing(mat, kwargs.pop('kNN', 20),
                                    kwargs.pop('iterations', 1), kwargs.pop('alpha', 0.9))
        
        elif op == 'logarithm':
            mat = np.log(mat + kwargs.pop('epsilon', 1)) / np.log(kwargs.pop('log_base', np.e))
        
        elif op == 'power':
            mat = np.power(mat, kwargs.pop('pow', 0.5))
        
    return mat


def OE_norm(mat):
    new_mat = np.zeros(mat.shape)
    averages = np.array([np.mean(mat[i:, :len(mat) - i]) for i in range(len(mat))])
    averages = np.where(averages == 0, 1, averages)
    for i in range(len(mat)):
        for j in range(len(mat)):
            d = abs(i - j)
            new_mat[i, j] = mat[i, j] / averages[d]
    return new_mat

def VC_norm(mat):
    sm = np.sum(mat, axis=0)
    sm = np.where(sm == 0, 1, sm)
    sm_v = np.tile(sm, (len(sm), 1))
    sm_c = sm_v.T
    new_mat = mat / sm_c / sm_v
    return new_mat

def VC_SQRT_norm(mat):
    sm = np.sum(mat, axis=0)
    sm = np.where(sm == 0, 1, sm)
    sm = np.sqrt(sm)
    sm_v = np.tile(sm, (len(sm), 1))
    sm_c = sm_v.T
    new_mat = mat / sm_c / sm_v
    return new_mat

def KR_norm(mat, maximum_error_rate=1e-4):
    bias = np.mean(mat) * maximum_error_rate
    # Remove all-zero rows and columns
    sm = np.sum(mat, axis=0)
    zeros = []
    for i in range(len(sm)):
        if sm[i] == 0:
            zeros.append(i)
    new_mat = np.delete(mat, zeros, axis=0)
    new_mat = np.delete(new_mat, zeros, axis=1)

    # Iteration
    x = np.random.random(size=len(new_mat))
    k = 0
    while True:
        # I forgot where I found this iteration formula
        # But it does work...
        # I'll check later...
        k += 1
        aa = np.diag(x).dot(new_mat) + np.diag(new_mat.dot(x))
        aa = np.linalg.inv(aa)
        bb = np.diag(x).dot(new_mat).dot(x) - np.ones(x.shape)
        delta = aa.dot(bb)
        new_x = x - delta
        max_error = np.max(np.abs(delta))
        # print(f'Iteration: {k}, Max Error: {max_error}')
        if max_error < bias:
            break
        else:
            x = new_x
    # Normalization
    dg = np.diag(new_x)
    new_mat = dg.dot(new_mat).dot(dg)
    # Put all-zero rows and columns back
    for zero in zeros:
        new_mat = np.insert(new_mat, zero, 0, axis=0)
        new_mat = np.insert(new_mat, zero, 0, axis=1)
    return new_mat

def convolution(mat, kernel_shape=3):
    conv = np.ones((kernel_shape, kernel_shape)) / (kernel_shape ** 2)
    mat = convolve2d(mat, conv, 'same')
    return mat

def random_walk(mat, random_walk_ratio=1.0,t=1):
    sm = np.sum(mat, axis=1)
    sm = np.where(sm == 0, 1, sm)
    sm = np.tile(sm, (len(mat), 1)).T
    walk = mat / sm
    for i in range(t):
        mat = random_walk_ratio * mat.dot(walk) + (1 - random_walk_ratio) * mat
    return mat

def reduce_sparsity(mat, sparsity_method='log', power=0.5):
    if sparsity_method == 'log':
        return np.log(mat + 1)
    elif sparsity_method == 'power':
        return np.power(mat, power)
    else:
        raise ValueError('Method {0} not supported while reducing sparsity.'.format(sparsity_method))

def network_enhancing(mat, kNN=20, iteration=1, alpha=0.9):
    argsort = np.argsort(-mat, axis=1)
    new_mat = np.zeros(mat.shape)
    for i in range(len(mat)):
        for j in range(kNN):
            pos = argsort[i, j]
            new_mat[i, pos] = mat[i, pos]
    sm = np.sum(new_mat, axis=1)
    sm = np.where(sm == 0, 1, sm)
    sm = np.tile(sm, (len(mat), 1)).T
    walk = new_mat / sm
    for k in range(iteration):
        if k == 0:
            new_mat = alpha * walk.T.dot(mat).dot(walk) + (1 - alpha) * mat
        else:
            new_mat = alpha * walk.T.dot(new_mat).dot(walk) + (1 - alpha) * new_mat
    return mat

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols
    
def oe_norm(matrix, expected = None):
    new_matrix = np.zeros_like(matrix)
    matrix =  np.array(matrix).astype('float32')
    for k in range(len(matrix)):
        rows, cols = kth_diag_indices(matrix, k)
        diag = np.diag(matrix,k)
        if expected is not None:
            expect = expected[k]
        else:
            expect = np.nanmean(diag)
        if expect == 0:
            new_matrix[rows, cols] = 0.0
        else:
            new_matrix[rows, cols] = diag / (expect+1e-15)
    new_matrix = new_matrix + new_matrix.T
    return new_matrix

def min_max_norm(mat):
    mat = (mat - mat.min())/(mat.max()-mat.min())
    return mat

def log_norm(mat):
    mat = np.log(mat + 1)
    return mat    

def get_matrix(cell,chrom='chr1'):
    matrix = cooler.Cooler(cell).matrix(balance=False).fetch(chrom)[:]
    matrix = np.log(oe_norm(matrix)+1)
    return matrix

