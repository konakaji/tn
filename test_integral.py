from unittest import TestCase
from computation import *
from core import *


class TestPathUtil(TestCase):
    def test_find_outside_pairs_1(self):
        lp1 = Plug(0, Direction.Left, 0)
        rp1 = Plug(0, Direction.Right, 0)
        lp2 = Plug(0, Direction.Left, 1)
        rp2 = Plug(0, Direction.Right, 1)
        lpd1 = Plug(0, Direction.Left, 2)
        rpd1 = Plug(0, Direction.Right, 2)
        lpd2 = Plug(0, Direction.Left, 3)
        rpd2 = Plug(0, Direction.Right, 3)

        ar = Plug(0, Direction.Right, 4)
        bl = Plug(0, Direction.Left, 5)
        cr = Plug(0, Direction.Right, 6)
        dl = Plug(0, Direction.Left, 7)
        edr = Plug(0, Direction.Right, 8)
        fdl = Plug(0, Direction.Left, 9)
        gdr = Plug(0, Direction.Right, 10)
        hdl = Plug(0, Direction.Left, 11)
        haar_pairs = [(lp1, rpd1), (lp2, rpd2), (lpd1, rp1), (lpd2, rp2)]
        edges = []
        edges.append(Edge(ar, lp1))
        edges.append(Edge(rp1, bl))
        edges.append(Edge(cr, lp2))
        edges.append(Edge(rp2, dl))
        edges.append(Edge(edr, lpd1))
        edges.append(Edge(rpd1, fdl))
        edges.append(Edge(gdr, lpd2))
        edges.append(Edge(rpd2, hdl))
        final_pairs, n_loop = PathUtil.find_outside_pairs(haar_pairs)
        self.assertEqual(0, n_loop)
        final_pairs_answer = [(ar, fdl), (edr, bl), (cr, hdl), (gdr, dl)]
        self.assertEquals(4, len(final_pairs))
        final_pairs = sorted(final_pairs, key=lambda v: v[0].node_id)
        final_pairs_answer = sorted(final_pairs_answer, key=lambda v: v[0].node_id)
        for i in range(len(final_pairs)):
            self.assertEquals(final_pairs[i][0], final_pairs_answer[i][0])
            self.assertEquals(final_pairs[i][1], final_pairs_answer[i][1])

    def test_find_outside_pairs_2(self):
        lp1 = Plug(0, Direction.Left, 0)
        rp1 = Plug(0, Direction.Right, 0)
        lp2 = Plug(0, Direction.Left, 1)
        rp2 = Plug(0, Direction.Right, 1)
        lpd1 = Plug(0, Direction.Left, 2)
        rpd1 = Plug(0, Direction.Right, 2)
        lpd2 = Plug(0, Direction.Left, 3)
        rpd2 = Plug(0, Direction.Right, 3)
        edges = []
        edges.append(Edge(rpd1, lp1))
        edges.append(Edge(rp1, lpd1))
        edges.append(Edge(rpd2, lp2))
        edges.append(Edge(rp2, lpd2))
        haar_pairs = [(lp1, rpd1), (lp2, rpd2), (lpd1, rp1), (lpd2, rp2)]
        final_pairs, n_loop = PathUtil.find_outside_pairs(haar_pairs)
        self.assertEqual(4, n_loop)
        self.assertEquals(0, len(final_pairs))

    def test_find_outside_pairs_3(self):
        lp1 = Plug(0, Direction.Left, 0)
        rp1 = Plug(0, Direction.Right, 0)
        lp2 = Plug(0, Direction.Left, 1)
        rp2 = Plug(0, Direction.Right, 1)
        lpd1 = Plug(0, Direction.Left, 2)
        rpd1 = Plug(0, Direction.Right, 2)
        lpd2 = Plug(0, Direction.Left, 3)
        rpd2 = Plug(0, Direction.Right, 3)
        edges = []
        edges.append(Edge(rpd2, lp1))
        edges.append(Edge(rp1, lpd2))
        edges.append(Edge(rpd1, lp2))
        edges.append(Edge(rp2, lpd1))
        haar_pairs = [(lp1, rpd1), (lp2, rpd2), (lpd1, rp1), (lpd2, rp2)]
        final_pairs, n_loop = PathUtil.find_outside_pairs(haar_pairs)
        self.assertEqual(2, n_loop)
        self.assertEquals(0, len(final_pairs))

    def test_find_outside_pairs_4(self):
        lp1 = Plug(0, Direction.Left, 0)
        rp1 = Plug(0, Direction.Right, 0)
        lp2 = Plug(0, Direction.Left, 1)
        rp2 = Plug(0, Direction.Right, 1)
        lpd1 = Plug(0, Direction.Left, 2)
        rpd1 = Plug(0, Direction.Right, 2)
        lpd2 = Plug(0, Direction.Left, 3)
        rpd2 = Plug(0, Direction.Right, 3)
        cr = Plug(0, Direction.Right, 6)
        dl = Plug(0, Direction.Left, 7)
        edr = Plug(0, Direction.Right, 8)
        fdl = Plug(0, Direction.Left, 9)
        edges = []
        edges.append(Edge(rpd2, lp1))
        edges.append(Edge(rp1, lpd2))
        edges.append(Edge(cr, lp2))
        edges.append(Edge(edr, lpd1))
        edges.append(Edge(rpd1, fdl))
        edges.append(Edge(rp2, dl))
        haar_pairs = [(lp1, rpd1), (lp2, rpd2), (lpd1, rp1), (lpd2, rp2)]
        final_pairs_answer = [(cr, fdl), (edr, dl)]
        final_pairs, n_loop = PathUtil.find_outside_pairs(haar_pairs)
        self.assertEqual(0, n_loop)
        self.assertEquals(2, len(final_pairs))
        final_pairs = sorted(final_pairs, key=lambda v: v[0].node_id)
        final_pairs_answer = sorted(final_pairs_answer, key=lambda v: v[0].node_id)
        for i in range(len(final_pairs)):
            self.assertEquals(final_pairs[i][0], final_pairs_answer[i][0])
            self.assertEquals(final_pairs[i][1], final_pairs_answer[i][1])

    def test_find_pairs(self):
        lp1 = Plug(0, Direction.Left, 0)
        rp1 = Plug(0, Direction.Right, 0)
        lp2 = Plug(0, Direction.Left, 1)
        rp2 = Plug(0, Direction.Right, 1)
        lpd1 = Plug(0, Direction.Left, 2)
        rpd1 = Plug(0, Direction.Right, 2)
        lpd2 = Plug(0, Direction.Left, 3)
        rpd2 = Plug(0, Direction.Right, 3)
        edges = []
        edges.append(Edge(rpd1, lp1))
        edges.append(Edge(rp1, lpd1))
        edges.append(Edge(rpd2, lp2))
        edges.append(Edge(rp2, lpd2))

        haar_pairs = [(lp1, rpd1), (lp2, rpd2), (lpd1, rp1), (lpd2, rp2)]
        self.assertEquals(4, len(PathUtil.find_pairs(haar_pairs)))
