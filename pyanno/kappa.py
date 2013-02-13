from .util import *

def agr(confusion_mat):
    """Return the agreement rate for the specified confusion matrix.
    
    Keyword arguments:
    confusion_mat -- square confusion matrix of catgorical responses
    """
    tot = mat_sum(confusion_mat)
    agr = 0
    k = len(confusion_mat)
    while k > 0:
        k -= 1
        agr += confusion_mat[k][k]
    return float(agr)/float(tot)


def s(confusion_mat):
    """Return the s statistic for the specified confusion matrix.

    Keyword arguments:
    confusion_mat -- square confusion matrix of catgorical responses
    """
    agr_ = agr(confusion_mat)
    e_agr = 1.0/float(len(confusion_mat))
    return chance_adj_agr(agr_,e_agr)


def pi(confusion_mat):
    """Return Scott's pi statistic for the specified confusion matrix.

    Keyword arguments:
    confusion_mat -- square confusion matrix of catgorical responses
    """
    agr_ = agr(confusion_mat)
    K = len(confusion_mat)
    Ks = range(K)
    theta = alloc_vec(K)
    for k1 in Ks:
        for k2 in Ks:
            theta[k1] += confusion_mat[k1][k2]
            theta[k2] += confusion_mat[k1][k2]
    prob_norm(theta)
    e_agr = 0.0
    for k in Ks:
        e_agr += theta[k]**2
    return chance_adj_agr(agr_,e_agr)


def kappa(confusion_mat):
    """Return Cohen's kappa statistic for the specified confusion matrix.

    Keyword arguments:
    confusion_mat -- square confusion matrix of catgorical responses
    """
    agr_ = agr(confusion_mat)
    K = len(confusion_mat)
    Ks = range(K)
    theta1 = alloc_vec(K)
    theta2 = alloc_vec(K)
    for k1 in Ks:
        for k2 in Ks:
            theta1[k1] += confusion_mat[k1][k2]
            theta2[k2] += confusion_mat[k1][k2]
    prob_norm(theta1)
    prob_norm(theta2)
    e_agr = 0.0
    for k in Ks:
        e_agr += theta1[k] * theta2[k]
    return chance_adj_agr(agr_,e_agr)


def chance_adj_agr(agr,expected_agr):
    """Return the chance-adjusted agreement given the specified agreement
    and expected agreement.  

    Defined by (agr - expected_agr)/(1.0 - expected_agr)

    Keyword arguments:
    agr -- agreement
    expected_agr -- expected agreement
    confusion_mat -- square confusion matrix of catgorical responses
    """
    return (agr - expected_agr)/(1.0 - expected_agr)
          

def K(item,anno,label):
    """Return the K agreement statistic for multiple annotators
    represented by the specified items, annotators, and labels.
    
    Keyword arguments:
    item -- array of item IDs, with item[n] being item in n-th annotation
    anno -- parallel array of annotator IDs, with anno[n] being the 
            annotator for the n-th annotation
    label -- parallel array of labels, with label[n] being the
             label assigned by annotator anno[n] to item item[n]
    """
    if len(item) != len(anno):
        raise ValueError("len(item) != len(anno)")
    if len(label) != len(anno):
        raise ValueError("len(label) != len(anno)")
    I = max(item)+1
    J = max(anno)+1
    K = max(label)+1
    N = len(item)
    theta = global_prevalence(item,label)
    agr_exp = 0.0
    for theta_n in theta:
        agr_exp += theta_n * theta_n
    anno_labels = alloc_vec(I)
    i = I
    while i > 0:
        i -= 1
        anno_labels[i] = []
    n = len(item)
    while n > 0:
        n -= 1
        anno_labels[item[n]].append(label[n])
    tot = 0
    agr = 0
    i = I
    while i > 0:
        i -= 1
        M = len(anno_labels[i])
        if M < 2:
            continue
        anno_labels[i].sort()
        tot += (M * (M-1)) / 2
        start = 0
        m = 1
        while m < M:
            if anno_labels[i][m] != anno_labels[i][start]:
                run = m - start
                agr += (run * (run-1)) / 2
                start = m
            m += 1
        run = m - start
        agr += (run * (run-1)) / 2
    return chance_adj_agr(float(agr)/float(tot),agr_exp)
            
            
def global_prevalence(item,label):
    """Return the global maximum likelihood estimate of prevalence for
    the specified array of labels.  

    The prevalence of a category is estimated to be proportional to
    the sum of adjusted counts for labels, where the adjusted count
    for label[n] is 1 divided by the number of labels for the item[n].

    The length of the returned array will be the number of labels, the
    values will be non-negative, and the values will sum to 1.0

    Keyword arguments:
    item -- array of item identifiers
    label -- parallel array of label assignments
    """
    if len(label) != len(item):
        raise ValueError("len(label) != len(item)")
    item_to_labs = {}
    for i in item:
        item_to_labs[i] = []
    N = len(item)
    n = 0
    while n < N:
        item_to_labs[item[n]].append(label[n])
        n += 1

    K = max(label)+1
    theta = alloc_vec(K)
    for labs in item_to_labs.values():
        increment = 1.0/float(len(labs))
        for k in labs:
            theta[k] += increment
    prob_norm(theta)
    return theta



