import sys
try:
    import pyanno.multinom
    import pyanno.kappa
    import pyanno.util
except ImportError, e:
    print e
    print ""
    print "Need to install pyanno package and dependencies."
    print "See instructions in Install.txt from pyanno distribution"
    raise SystemExit(1)

if len(sys.argv) != 2:
    print "Require command-line argument for data file."
    print "Data is distributed in $PYANNO_HOME/examples/rzhetsky_2009/data"
    raise SystemExit(2)

print "READING DATA"
file_name = sys.argv[1]
print "     input file=",file_name
file = open(file_name)
item = []
label = []
anno = []
item_id = 0

while True:
    line = file.readline()
    if not line:
        break
    fields = line.split(',')
    if len(fields) != 13:
        continue
    for j in range(13):
        field_int = int(fields[j])
        if field_int > 0:
            item.append(item_id)
            label.append(field_int - 1)
            anno.append(j)
    item_id += 1
file.close()

I = item_id
J = max(anno)+1
K = max(label)+1
N = len(item)

Is = range(I)
Js = range(J)
Ks = range(K)
Ns = range(N)

print "SIZES"
print "    I=",I
print "    J=",J
print "    K=",K
print "    N=",N

print "VOTED PREVALENCE"
print "    prev=",pyanno.kappa.global_prevalence(item,label)

print "PRIORS"
beta = pyanno.util.alloc_vec(K,1.5)  
alpha = pyanno.util.alloc_mat(K,K,1.5)
for k1 in Ks:
    for k2 in Ks:
        if (k1 == k2):
            alpha[k1][k2] = 8.0
            if k2 - 1 >= 0:
                alpha[k1][k2-1] = 2.0
            if k2 + 1 < K:
                alpha[k1][k2+1] = 2.0
print "     beta=",beta
for k in Ks:
    print "     alpha[",k,"]=",alpha[k]

print "RUNNING EM"
(diff,ll,lp,prev_map,cat_map,accuracy_map) \
    = pyanno.multinom.map(item,anno,label,alpha,beta,init_acc=0.6,epsilon=0.000001)

print "CONVERGENCE ll[final] - ll[final-10]=",diff

print "PREVALENCE ESTIMATES"
print "{0:>2s}, {1:>5s}".format("k","MAP")
for k in Ks:
    print "{0:2d}, {1:5.3f}".format(k,prev_map[k])

print "ACCURACY ESTIMATES"
print "{0:>3s},{1:>2s},{2:>2s}, {3:>5s}".format("j","k1","k2","map")
for j in Js:
    for k1 in Ks:
        for k2 in Ks:
            print "{0:3d},{1:2d},{2:2d}, {3:5.3f}".\
                format(j,k1,k2,accuracy_map[j][k1][k2])

print "CATEGORY ESTIMATES"
for i in Is:
    print "{0:5d}".format(i),
    for k in Ks:
        print " {0:2d},{1:5.3f}".format(k,cat_map[i][k]),
    print ""

