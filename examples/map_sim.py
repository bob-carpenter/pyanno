try:
    import pyanno.multinom
    import pyanno.util
except ImportError, e:
    print e
    print ""
    print "Need to install pyanno package and dependencies."
    print "See instructions in Install.txt from pyanno distribution"
    raise SystemExit(1)
    
# Simulated Sizes
I = 200
J = 5
K = 4
N = I*J 

Is = range(I)
Js = range(J)
Ks = range(K)
Ns = range(N)

print "SIZES"
print "I=",I
print "J=",J
print "K=",K
print "N=",N

print "SIMULATING ORDINAL CODING DATA SET"
(prev,cat,accuracy,item,anno,label) = pyanno.multinom.sim_ordinal(I,J,K)

print "CALCULATING SAMPLE PARAMETERS"
prev_sample = pyanno.multinom.alloc_vec(K)
for k in cat:
    prev_sample[k] += 1
pyanno.multinom.prob_norm(prev_sample)

accuracy_sample = pyanno.multinom.alloc_tens(J,K,K)
for n in Ns:
    accuracy_sample[anno[n]][cat[item[n]]][label[n]] += 1
for j in Js:
    for k in Ks:
        pyanno.multinom.prob_norm(accuracy_sample[j][k])

alpha = pyanno.util.alloc_mat(K,K,1.0) # add 1 throughout
for k in Ks:
    alpha[k][k] = 4
    if k-1 >= 0:
        alpha[k][k-1] = 2
    if k+1 < K:
        alpha[k][k+1] = 2
beta = pyanno.util.alloc_vec(K,1.0) # add 1 smoothing

print "PRIORS"
print "    alpha=",alpha
print "    beta=",beta

print "RUNNING EM"
epsilon = 0.00001
init_acc = 0.6
(diff,ll,lp,prev_map,
 cat_map,accuracy_map) = pyanno.multinom.map(item,anno,label,
                                             alpha,beta,init_acc,epsilon)


print "CONVERGENCE ll[final] - ll[final-10]=",diff

print "PREVALENCE ESTIMATES"
print "{0:>2s}, {1:>5s}, {2:>5s}, {3:>5s}, {4:>6s}, {5:>6s}".\
        format("k","sim","samp","MAP","d.sim","d.samp")
for k in Ks:
    print "{0:2d}, {1:5.3f}, {2:5.3f}, {3:5.3f}, {4:+5.3f}, {5:+5.3f}".\
            format(k,prev[k],prev_sample[k],prev_map[k],
                   prev_map[k]-prev[k],
                   prev_map[k]-prev_sample[k])

print "ACCURACY ESTIMATES"
print "{0:>3s},{1:>2s},{2:>2s}, {3:>5s}, {4:>5s}, {5:>5s}, {6:>6s}, {7:>6s}".\
        format("j","k1","k2","sim","samp","map","d.sim","d.samp")
for j in Js:
    for k1 in Ks:
        for k2 in Ks:
            print "{0:3d},{1:2d},{2:2d}, {3:5.3f}, {4:5.3f}, {5:5.3f}, {6:+5.3f}, {7:+5.3f}".\
                    format(j,k1,k2,accuracy[j][k1][k2],
                           accuracy_sample[j][k1][k2],
                           accuracy_map[j][k1][k2],
                           accuracy_map[j][k1][k2]-accuracy[j][k1][k2],
                           accuracy_map[j][k1][k2]-accuracy_sample[j][k1][k2])

print "CATEGORY ESTIMATES"
for i in Is:
    print "{0:5d}".format(i),
    for k in Ks:
        match = "*" if cat[i] == k else " "
        print " {0:2d}:{1:1}{2:5.3f}".format(k,match,cat_map[i][k]),
    print ""
