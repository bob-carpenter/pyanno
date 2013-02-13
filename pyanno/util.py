def alloc_vec(N,x=0.0):
    result = []
    n = 0
    while n < N:
        result.append(x)
        n += 1
    return result


def alloc_mat(M,N,x=0.0):
    result = []
    m = 0
    while m < M:
        result.append(alloc_vec(N,x))
        m += 1
    return result
                   
    
def alloc_tens(M,N,J,x=0.0):
    result = []
    for m in range(M):
        result.append(alloc_mat(N,J,x))
    return result


def alloc_tens4(M,N,J,K,x=0.0):
    result = []
    for m in range(M):
        result.append(alloc_tens(N,J,K,x))
    return result

            
def fill_vec(xs,y):
    i = 0
    while i < len(xs):
        xs[i] = y
        i += 1


def fill_mat(xs,y):
    i = 0
    while i < len(xs):
        fill_vec(xs[i],y)
        i += 1


def fill_tens(xs,y):
    i = 0
    while i < len(xs):
        fill_mat(xs[i],y)
        i += 1


def vec_copy(x,y):
    n = len(x)
    while (n > 0):
        n -= 1
        y[n] = x[n]


def vec_sum(x):
    sum = 0
    for x_i in x:
        sum += x_i
    return sum


def mat_sum(x):
    sum = 0
    for x_i in x:
        sum += vec_sum(x_i)
    return sum


def prob_norm(theta):
    Z = sum(theta)
    if Z <= 0.0:
        fill_vec(theta,1.0/float(len(theta)))
        return
    n = len(theta) - 1
    while n >= 0:
        theta[n] /= Z
        n -= 1


def warn_missing_vals(varname,xs):
    missing = set(xs) - set(range(max(xs)+1))
    if len(missing) > 0:
        print "Missing values in ",varname,"=",missing



        
    
        
    
