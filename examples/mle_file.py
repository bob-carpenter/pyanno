import sys
try:
    import pyanno.multinom
    import pyanno.util
    import pyanno.kappa
except ImportError, e:
    print e
    print ""
    print "Need to install pyanno package and dependencies."
    print "See instructions in Install.txt from pyanno distribution"
    raise SystemExit(1)

#if len(sys.argv) != 2:
#    print "Require name of data file as command-line argument."

# return (number-of-symbols, id-list, id-to-item-list, item-to-id-dict)
def to_sym_tab(xsyms):
    sym2id = {}
    id2sym = []
    xs = []
    next_id = 0
    for xsym in xsyms:
        if xsym not in sym2id:
            sym2id[xsym] = next_id
            id2sym.append(xsym)
            next_id += 1
        xs.append(sym2id[xsym])
    return (next_id,xs,id2sym,sym2id)
        
print "READING DATA"
file_name = sys.argv[1]
print "     input file=",file_name
item_sym = []
coder_sym = []
label_sym = []
file = open(file_name)
line_num = 1
while True:
    line = file.readline()
    if not line:
        break
    fields = line.split(',')
    if len(fields) != 3:
        print "Ignoring ill-formed line number ",line_num
        print "     line=|",line,"|"
        continue
    item_sym.append(fields[0])
    coder_sym.append(fields[1])
    label_sym.append(fields[2])
    line_num += 1
file.close()
N = len(item_sym)
(I,item,item_id2sym,item_sym2id) = to_sym_tab(item_sym)
(J,coder,coder_id2sym,coder_sym2id) = to_sym_tab(coder_sym)
(K,label,label_id2sym,label_sym2id) = to_sym_tab(label_sym)

print "VOTED PREVALENCE"
print "    prev=",pyanno.kappa.global_prevalence(item,label)

print "RUNNING EM"
(diff,ll,prev_mle,cat_mle,accuracy_mle) \
    = pyanno.multinom.mle(item,coder,label,init_acc=0.6,epsilon=0.000001)

print "CONVERGENCE ll[final] - ll[final-10]=",diff

print "PREVALENCE ESTIMATES"
print "{0:>2s}, {1:>5s}".format("k","MLE")
for k in range(K):
    print "{0:2d}, {1:5.3f}".format(k,prev_mle[k])

print "ACCURACY ESTIMATES"
print "{0:>3s},{1:>2s},{2:>2s}, {3:>5s}".format("j","k1","k2","mle")
for j in range(J):
    for k1 in range(K):
        for k2 in range(K):
            print "{0:3d},{1:2d},{2:2d}, {3:5.3f}".\
                format(j,k1,k2,accuracy_mle[j][k1][k2])

print "CATEGORY ESTIMATES"
for i in range(I):
    print "{0:5d}".format(i),
    for k in range(K):
        print " {0:2d},{1:5.3f}".format(k,cat_mle[i][k]),
#        print " {0:8s},{1:5.3f}".format(item_id2sym[k],cat_mle[i][k]),
    print ""

print "ITEMS"
c = 0
for xyz in item_id2sym:
    print c, xyz
    c+=1
