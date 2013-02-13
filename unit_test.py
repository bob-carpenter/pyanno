from pyanno.kappa import *
from pyanno.multinom import *
from pyanno.util import *

import unittest

class TestUtil(unittest.TestCase):
    def test_vec_copy(self):
        c = []
        d = []
        vec_copy(c,d)
        self.assertEquals(c,d)
        a = [1,2,3]
        b = [4,5,6]
        vec_copy(a,b)
        self.assertEquals(a,b)
    def test_fill_vec(self):
        a = []
        fill_vec(a,0)
        self.assertEquals([],a)
        b = [1]
        fill_vec(b,3)
        self.assertEquals([3],b)
        c = [1,2,3]
        fill_vec(c,5)
        self.assertEquals([5,5,5],c)
    def test_fill_mat(self):
        a = [[]]
        fill_mat(a,0)
        self.assertEquals([[]],a)
        b = [[1,2,3]]
        fill_mat(b,5)
        self.assertEquals([[5,5,5]],b)
        c = [[1,2],[3,4]]
        fill_mat(c,6)
        self.assertEquals([[6,6],[6,6]],c)
    def test_fill_tens(self):
        a = [[[]]]
        fill_tens(a,0)
        self.assertEquals([[[]]],a)
        b = [[[1,2],[3,4]],[[5,6],[7,8]]]
        fill_tens(b,10)
        self.assertEquals([[[10,10],[10,10]],[[10,10],[10,10]]],b)
    def test_alloc_vec(self):
        a = alloc_vec(0)
        self.assertEquals([],a)
        b = alloc_vec(1)
        self.assertEquals([0.0],b)
        c = alloc_vec(2,9.0)
        self.assertEquals([9.0,9.0],c)
    def test_alloc_mat(self):
        a = alloc_mat(1,0)
        self.assertEquals([[]],a)
        b = alloc_mat(1,1,3.0)
        self.assertEquals([[3.0]],b)
        c = alloc_mat(2,3,5.0)
        self.assertEquals([[5.0,5.0,5.0],[5.0,5.0,5.0]],c)
    def test_alloc_tens(self):
        a = alloc_tens(1,1,0)
        self.assertEquals([[[]]],a)
        b = alloc_tens(1,1,1)
        self.assertEquals([[[0.0]]],b)
        c = alloc_tens(2,3,4,9)
        self.assertEquals([[[9,9,9,9],[9,9,9,9],[9,9,9,9]],
                           [[9,9,9,9],[9,9,9,9],[9,9,9,9]]],c)
    def test_alloc_tens4(self):
        a = alloc_tens4(1,1,1,0)
        self.assertEquals([[[[]]]],a)
        b = alloc_tens4(1,1,1,1)
        self.assertEquals([[[[0.0]]]],b)
        c = alloc_tens4(1,2,3,4,9)
        self.assertEquals([[[[9,9,9,9],[9,9,9,9],[9,9,9,9]],
                            [[9,9,9,9],[9,9,9,9],[9,9,9,9]]]],c)
    def test_vec_sum(self):
        self.assertEquals(0,vec_sum([]))
        self.assertEquals(1,vec_sum([1]))
        self.assertEquals(3,vec_sum([1,2]))
    def test_mat_sum(self):
        self.assertEquals(0,mat_sum([[]]))
        self.assertEquals(1,mat_sum([[1]]))
        self.assertEquals(3,mat_sum([[1,2]]))
        self.assertEquals(21,mat_sum([[1,2,3],[4,5,6]]))
    def test_prob_norm(self):
        theta = [0.2]
        prob_norm(theta)
        self.assert_prob_normed(theta)


    def assert_prob_normed(self,theta):
        self.assert_(len(theta) > 0)
        for theta_i in theta:
            self.assert_(theta_i >= 0.0)
            self.assert_(theta_i <= 1.0)
        self.assertAlmostEqual(1.0,vec_sum(theta),3)

class TestKappa(unittest.TestCase):
    def test_agr(self):
        cm1 = [[1]]
        self.assertEquals(1.0,agr(cm1))
        cm2 = [[41,3],[9,47]]
        self.assertAlmostEqual((41.0 + 47.0)/100.0,agr(cm2))
        cm3 = [[44,6],[6,44]]
        self.assertAlmostEqual(0.88,agr(cm3))
    def test_s(self):
        cm1 = [[44,6],[6,44]]
        self.assertAlmostEqual((0.88-0.5)/(1.0-0.5),s(cm1))
        cm2 = [[44,6,0,0],[6,44,0,0],[0,0,0,0],[0,0,0,0]]
        self.assertAlmostEqual(0.84,s(cm2))
        self.assertAlmostEqual(0.88,agr(cm2))
    def test_pi(self):
        cm1 = [[44,6,0],[6,44,0],[0,0,0]]
        self.assertAlmostEqual(0.88,agr(cm1))
        self.assertAlmostEqual(0.82,s(cm1))
        self.assertAlmostEqual(0.76,pi(cm1))
        cm2 = [[77,1,2],[1,6,3],[2,3,5]]
        self.assertAlmostEqual(0.88,agr(cm2))
        self.assertAlmostEqual(0.82,s(cm2))
        self.assertAlmostEqual((0.88-0.66)/(1.0-0.66),pi(cm2))
        cm3 = [[990,5],[5,0]]
        self.assertAlmostEqual(0.99,agr(cm3))
        self.assertAlmostEqual(0.98,s(cm3))
        self.assertAlmostEqual((0.99-0.99005)/(1.0-0.99005),pi(cm3))
    def test_kappa(self):
        cm1 = [[8,1],
               [4,2]]
        ce = ((12.0*9.0) + (6.0*3.0))/(15.0*15.0)
        ag = 10.0/15.0
        k = (ag - ce)/(1.0 - ce)
        self.assertAlmostEqual(k,kappa(cm1))
        cm2 = [[7,5,0],
               [1,9,2],
               [0,8,4]]
        ce2 = (8.0*12.0 + 22.0*12.0 + 6.0*12.0)/(36.0*36.0)
        ag2 = 20.0/36.0
        k2 = (ag2 - ce2)/(1.0 - ce2)
        self.assertAlmostEqual(k2,kappa(cm2))
    def test_global_prevalence(self):
        item = [0,0]
        label = [0,1]
        theta_hat = global_prevalence(item,label)
        assertAlmostEqualLists(self,[0.5,0.5],theta_hat)
        item2 = [0,0,0, 1,1]
        label2 = [0,0,1, 1,2]
        theta_hat2 = global_prevalence(item2,label2)
        assertAlmostEqualLists(self,[4.0/12.0, 5.0/12.0, 3.0/12.0],theta_hat2)
    def test_K(self):
        item = [0,0,0]
        anno = [0,1,2]
        lab = [0,0,1]
        ag = 1.0/3.0
        ag_exp = (1.0*1.0)/(3.0*3.0) + (2.0*2.0)/(3.0*3.0)
        Ke = (ag - ag_exp)/(1.0 - ag_exp)
        self.assertAlmostEqual(Ke,K(item,anno,lab))
    def test_K_2(self):
        item = [0,0,0, 1,1, 2]
        anno = [0,1,2, 0,1, 0]
        lab =  [0,0,1, 1,1, 0]
        ag = 2.0/4.0
        ag_exp = (5.0/9.0)**2 + (4.0/9.0)**2
        Ke = (ag - ag_exp)/(1.0 - ag_exp)
        self.assertAlmostEqual(Ke,K(item,anno,lab))
    def test_K_3(self):
        item = [0,0,0, 1,1, 2,2,2,2,2,2,2]
        anno = [0,1,2, 0,1, 0,1,2,3,4,5,6]
        lab =  [0,0,1, 1,1, 0,1,2,1,0,1,2]
        ag = (1.0 + 1.0 + 1.0 + 3.0 + 1.0)/(3.0 + 1.0 + 21.0)
        ag_exp = (20.0/63.0)**2 + (37.0/63.0)**2 + (6.0/63.0)**2
        Ke = (ag - ag_exp)/(1.0 - ag_exp)
        self.assertAlmostEqual(Ke,K(item,anno,lab))
    def test_chance_adj_agr(self):
        self.assertAlmostEqual((0.9-0.5)/(1.0-0.5),chance_adj_agr(0.9,0.5))

class TestMultinom(unittest.TestCase):
    def testbase2(self):
        self.assertEquals(2,2)

def assertAlmostEqualLists(testcase,x,y):
    testcase.assertEquals(len(x),len(y))
    n = 0
    while n < len(x):
        testcase.assertAlmostEqual(x[n],y[n])
        n += 1

if __name__ == '__main__':
    unittest.main()
