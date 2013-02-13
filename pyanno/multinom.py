# Library Module
import math
import numpy
import pymc
from .util import *

def mle(item,
        anno,
        label,
        init_acc=0.6,
        epsilon=0.0001,
        max_epochs=500):
    """Returns a tuple (diff,ll,prev,cat,acc) containing the final
    difference in log likelihood, estimates for the data's
    log likelihood, the prevalence of categories, the categories of
    items, and the accuracies of each annotator (see below for more
    detail, including dimensions and types).

    The model is defined in the Models.txt file in the top level of
    the distribution.

    The data is provided in the form of three parallel arrays of
    length N.  For each index n, the parallel arrays has the content

    item[n] is the ID of the item being annotated, 
    anno[n] the ID of the annotator doing the annotation, and
    label[n] the ID of the category the annotator assigned the item.
    
    The rest of this description assumes the following constants:

    I = max(item) + 1 is the number of items, numbered 0 to I-1.
    J = max(anno) + 1 is the number of annotators, numbered 0 to J-1.
    K = max(label) + 1 is the number of categories, numbered 0 to K-1.

    Keyword arguments:
    item -- int[N], 0 <= item[n] < I
    list of item identifiers

    anno -- int[N], 0 <= anno[n] < J
    list of annotator identifiers

    label -- int[N], 0 <= label[n] < K
    list of labels by annotators for items

    init_acc --  float[J][K][K] >= 0.0, SUM_k2 init_acc[j][k][k2] = 1.0
    initial accuracy estimate for annotators for EM (default 0.6)

    epsilon -- float >= 0.0
    minimum log likelihood increase to continue (default 0.0001)

    max_epochs -- int > 0
    maximum numbe of epochs (default 500)

    Return tuple: (diff, ll, mle, cat, acc)
    diff --   float >= 0.0
    increase in log likelihood over past 10 epochs
    
    ll -- float <= 0.0
    log likelihood for estimates returned

    prev -- float[K] >= 0, SUM_k prev_mle[k] = 1.0
    maximum likelihood estimate of category prevalence

    cat -- float[I][K] >=0,  SUM_k cat_mle[i][k] = 1.0
    maximum likelihood estimate of categories

    acc -- float[J][K][K] >=0, SUM_k2 acc_mle[i][k][k2] = 1.0
    maximum likelihood estimate of annotator accuracies 
    """
    if epsilon < 0.0:
        raise ValueError("epislon < 0.0")
    if max_epochs < 0:
        raise ValueError("max_epochs < 0")

    log_likelihood_curve = []
    epoch = 0
    diff = float('inf')
    for (ll,prev_mle,cat_mle,accuracy_mle) in mle_em(item,anno,label,init_acc):
        print "  epoch={0:6d}  log lik={1:+10.4f}   diff={2:10.4f}".\
                format(epoch,ll,diff)
        log_likelihood_curve.append(ll)
        if epoch > max_epochs:
            break
        if len(log_likelihood_curve) > 10:
            diff = (ll - log_likelihood_curve[epoch-10])/10.0
            if abs(diff) < epsilon:
                break
        epoch += 1
    return (diff,ll,prev_mle,cat_mle,accuracy_mle)


def mle_em(item,
           anno,
           label,
           init_accuracy=0.6):
    """Returns a generator over per-epoch maximum-likelihood estimates,
    with generated items consisting of (ll,prev,cat,acc) tuples.

    See the documentation of pyanno.multinom.mle() for more information
    on the inputs.
    
    Keyword arguments:
    item -- int[N], 0 <= item[n] < I
    list of item identifiers

    anno -- int[N], 0 <= anno[n] < J
    list of annotator identifiers

    label -- int[N], 0 <= label[n] < K
    list of labels by annotators for items

    init_accuracy -- float, 0 <= init_accuracy <= 1
    initial accuracy of annotators to seed EM

    Generated tuples (ll,prev,cat,acc):
    ll -- log likelihood of data given parameters in
    current epoch
    prev -- estimate of prevalence in current epoch
    cat -- estimate of category distribution for items
    in current epoch
    acc -- estimate of annotator accuracies in current
    epoch
    """
    I = max(item) + 1
    J = max(anno) + 1
    K = max(label) + 1
    N = len(item)

    Is = range(I)
    Js = range(J)
    Ks = range(K)
    Ns = range(N)

    if len(anno) != N:
        raise ValueError("len(item) != len(anno)")
    if len(label) != N:
        raise ValueError("len(item) != len(label)")
    if init_accuracy < 0.0 or init_accuracy > 1.0:
        raise ValueError("init_accuracy not in [0,1]")
    for n in Ns:
        if item[n] < 0:
            raise ValueError("item[n] < 0")
        if anno[n] < 0:
            raise ValueError("anno[n] < 0")
        if label[n] < 0:
            raise ValueError("label[n] < 0")

    warn_missing_vals("item",item)
    warn_missing_vals("anno",anno)
    warn_missing_vals("label",label)

    # initialize params
    prevalence = alloc_vec(K,1.0/K)
    category = alloc_mat(I,K,1.0/K)
    accuracy = alloc_tens(J,K,K,(1.0 - init_accuracy)/(K-1.0))
    for j in Js:
        for k in Ks:
            accuracy[j][k][k] = init_accuracy

    
    while True:

        # E: p(cat[i]|...) 
        for i in Is:
            vec_copy(prevalence,category[i])
        for n in Ns:
            for k in Ks:
                category[item[n]][k] *= accuracy[anno[n]][k][label[n]]

        # log likelihood here to reuse intermediate category calc
        log_likelihood = 0.0
        for i in Is:
            likelihood_i = 0.0
            for k in Ks:
                likelihood_i += category[i][k]
            log_likelihood_i = math.log(likelihood_i)
            log_likelihood += log_likelihood_i

        for i in Is:
            prob_norm(category[i])


        # return here with E[cat|prev,acc] and LL(prev,acc;y)
        yield (log_likelihood,prevalence,category,accuracy)

        # M: prevalence* + accuracy*
        fill_vec(prevalence,0.0)
        for i in Is:
            for k in Ks:
                prevalence[k] += category[i][k]
        prob_norm(prevalence)

        fill_tens(accuracy,0.0)
        for n in Ns:
            for k in Ks:
                accuracy[anno[n]][k][label[n]] += category[item[n]][k]
        for j in Js:
            for k in Ks:
                prob_norm(accuracy[j][k])


def map(item,
        anno,
        label,
        alpha=None,
        beta=None,
        init_acc=0.6,
        epsilon=0.001,
        max_epochs=1000):
    """Returns maximum a posteriori (MAP) estimate tuple
    (diff,ll,lp,prev,cat,acc) consisting of final difference,
    log likelihood, log prior, category prevalence estimate,
    item category estimates, and annotator accuracy estimates.

    The model is defined in the Models.txt file in the top level of
    the distribution.
    
    See the documentation for pyanno.multinom.mle() in this module for
    a description of all but the following inputs:
    
    alpha -- matrix of priors for annotator accuracy by category
    (default uniform with all 1.0 values)

    beta -- array of priors for prevalence (default
    uniform with all 1.0 values)

    The outputs are the same as for pyanno.multinom.mle(), with the
    following addition:

    lp -- log of the prior p(acc|alpha) * p(prev|beta)
    """
    if epsilon < 0.0:
        raise ValueError("epislon < 0.0")
    if max_epochs < 0:
        raise ValueError("max_epochs < 0")

    llp_curve = []
    epoch = 0
    diff = float('inf')
    for (lp,ll,prev_mle,cat_mle,accuracy_mle) in map_em(item,anno,label,
                                                        alpha,beta,init_acc):
        print "  epoch={0:6d}  log lik={1:+10.4f}  log prior={2:+10.4f}  llp={3:+10.4f}   diff={4:10.4f}".\
                format(epoch,ll,lp,ll+lp,diff)
        llp_curve.append(ll+lp)
        if epoch > max_epochs:
            break
        if len(llp_curve) > 10:
            diff = (llp_curve[epoch] - llp_curve[epoch-10])/10.0
            if abs(diff) < epsilon:
                break
        epoch += 1
    return (diff,ll,lp,prev_mle,cat_mle,accuracy_mle)


def map_em(item,
           anno,
           label,
           alpha=None,
           beta=None,
           init_accuracy=0.6):
    """
    Return a generator of tuples (lp,ll,prev,cat,acc), one per epoch,
    consisting of log prior, log likelihood, prevalence, category, and
    accuracy estimates.  

    The MAP model is defined in Models.txt in the top-level directory.

    The arguments are the same as for pyanno.multinom.map(), and
    the tuples are just per-epoch versions of the same.

    This generator is to map() as mle_em() is to mle().  See the
    documentation for those functions in this module for more info.
    """
    I = max(item)+1
    J = max(anno)+1
    K = max(label)+1
    N = len(item)

    Is = range(I)
    Js = range(J)
    Ks = range(K)
    Ns = range(N)

    if alpha == None:
        alpha = alloc_mat(K,K,1.0)
    if beta == None:
        beta = alloc_vec(K,1.0)

    for k in Ks:
        if beta[k] < 1.0:
            raise ValueError("beta[k] < 1")
    for k1 in Ks:
        for k2 in Ks:
            if alpha[k1][k2] < 1.0:
                raise ValueError("alpha[k1][k2] < 1")

    alpha_prior_count = alloc_mat(K,K)
    for k1 in Ks:
        for k2 in Ks:
            alpha_prior_count[k1][k2] = alpha[k1][k2] - 1.0
    beta_prior_count = alloc_vec(K)
    for k in Ks:
        beta_prior_count[k] = beta[k] - 1.0

    beta_array = numpy.array(beta)
    alpha_array = []
    for k in Ks:
        alpha_array.append(numpy.array(alpha[k]))


    if len(anno) != N:
        raise ValueError("len(item) != len(anno)")
    if len(label) != N:
        raise ValueError("len(item) != len(label)")
    if init_accuracy < 0.0 or init_accuracy > 1.0:
        raise ValueError("init_accuracy not in [0,1]")
    for n in Ns:
        if item[n] < 0:
            raise ValueError("item[n] < 0")
        if anno[n] < 0:
            raise ValueError("anno[n] < 0")
        if label[n] < 0:
            raise ValueError("label[n] < 0")
    if len(alpha) != K:
        raise ValueError("len(alpha) != K")
    for k in Ks:
        if len(alpha[k]) != K:
            raise ValueError("len(alpha[k]) != K")
    if len(beta) != K:
        raise ValueError("len(beta) != K")

    warn_missing_vals("item",item)
    warn_missing_vals("anno",anno)
    warn_missing_vals("label",label)

    # initialize params
    prevalence = alloc_vec(K)
    vec_copy(beta_prior_count,prevalence)
    prob_norm(prevalence)

    category = alloc_mat(I,K)
    for i in Is:
        vec_copy(beta_prior_count,category[i])
        prob_norm(category[i])

    accuracy = alloc_tens(J,K,K,(1.0 - init_accuracy)/(K-1.0))
    for j in Js:
        for k in Ks:
            accuracy[j][k][k] = init_accuracy
    
    while True:

        # E: p(cat[i]|...) 
        for i in Is:
            vec_copy(prevalence,category[i])
        for n in Ns:
            for k in Ks:
                category[item[n]][k] *= accuracy[anno[n]][k][label[n]]

        # need log p(prev|beta) + SUM_k log p(acc[k]|alpha[k])
        # log likelihood here to reuse intermediate category calc
        log_likelihood = 0.0
        for i in Is:
            likelihood_i = 0.0
            for k in Ks:
                likelihood_i += category[i][k]
            if likelihood_i < 0.0:
                print "likelihood_i=",likelihood_i, "cat[i]=",category[i]
            log_likelihood_i = math.log(likelihood_i)
            log_likelihood += log_likelihood_i

        log_prior = 0.0
        prevalence_a = numpy.array(prevalence[0:(K-1)])
        log_prior += dir_ll(prevalence_a,beta_array)
        for j in Js:
            for k in Ks:
                acc_j_k_a = numpy.array(accuracy[j][k][0:(K-1)])
                log_prior += dir_ll(acc_j_k_a,alpha_array[k])
        if math.isnan(log_prior) or math.isinf(log_prior):
            log_prior = 0.0
        

        for i in Is:
            prob_norm(category[i])

        # return here with E[cat|prev,acc] and LL(prev,acc;y)
        yield (log_prior,log_likelihood,prevalence,category,accuracy)

        # M: prevalence* + accuracy*
        vec_copy(beta_prior_count,prevalence)
        for i in Is:
            for k in Ks:
                prevalence[k] += category[i][k]
        prob_norm(prevalence)

        for j in Js:
            for k in Ks:
                vec_copy(alpha_prior_count[k],accuracy[j][k])
        for n in Ns:
            for k in Ks:
                accuracy[anno[n]][k][label[n]] += category[item[n]][k]
        for j in Js:
            for k in Ks:
                prob_norm(accuracy[j][k])


# defined to prevent underflows resulting from theta[k] = 0.0, causing nans:
# >>> theta = numpy.array([ 0.75300156,  0.24474181,  0.00225663])
# >>> alpha = numpy.array([4.0,2.0,1.0,1.0])
# >>> pymc.dirichlet_like(theta,alpha)
# nan
def dir_ll(theta,alpha):
    delta = 0.0000000001
    while True:
        for k in range(len(theta)):
            # subtract delta, but with delta min
            theta[k] = max(delta, theta[k] - delta) 
        ll = pymc.dirichlet_like(theta,alpha)
        if not math.isnan(ll) and not math.isinf(ll):
            return ll

def sim_ordinal(I,J,K,alpha=None,beta=None):

    # test input params here

    Is = range(I)
    Js = range(J)
    Ks = range(K)
    N = I*J
    Ns = range(N)
        
    if alpha == None:
        alpha = alloc_mat(K,K)
        for k1 in Ks:
            for k2 in Ks:
                alpha[k1][k2] = max(1,(K + (0.5 if k1 == k2 else 0) \
                                       - abs(k1 - k2))**4)
        
    if beta == None:
        beta = alloc_vec(K,2.0)

    # simulated params
    beta = alloc_vec(K,2.0)

    prevalence = pymc.rdirichlet(beta).tolist()
    prevalence.append(1.0-sum(prevalence)) # complete
    category = []
    for i in Is:
        category.append(pymc.rcategorical(prevalence).tolist())

    accuracy = alloc_tens(J,K,K)
    for j in Js:
        for k in Ks:
            accuracy[j][k] = pymc.rdirichlet(alpha[k]).tolist()
            accuracy[j][k].append(1.0-sum(accuracy[j][k]))

    # simulated data
    item = []
    anno = []
    label = []
    for i in Is:
        for j in Js:
            item.append(i)
            anno.append(j)
            label.append(pymc.rcategorical(accuracy[j][category[i]]).tolist())
    N = len(item)

    return (prevalence,category,accuracy,item,anno,label)



