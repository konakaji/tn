from core import *
from abc import ABC, abstractmethod


class TNComputation(ABC):
    @abstractmethod
    def compute(self, networks: TensorNetworks) -> TensorNetworks:
        return networks


class Multiply(TNComputation):
    def __init__(self, type_left: Type, type_right: Type, left_dagger, right_dagger,
                 type_result: Type, result_dagger):
        self.type_left = type_left
        self.type_right = type_right
        self.left_dagger = left_dagger
        self.right_dagger = right_dagger
        self.type_result = type_result
        self.result_dagger = result_dagger

    def compute(self, networks: TensorNetworks):
        for network in networks.networks:
            self.do_compute(network)
        return networks

    def do_compute(self, network: TensorNetwork):
        pairs = []
        for n in network.nodes():
            n: Gate = n
            if n.type != self.type_left or n.dagger != self.left_dagger:
                continue
            n2 = None
            for p in n.get_right_plugs():
                node = network.node_map[p.edge.right_plug.node_id]
                if n2 is not None and n2 != node:
                    n2 = None
                    break
                if node.type != self.type_right or node.dagger != self.right_dagger:
                    break
                n2 = node
            if n2 is not None:
                pairs.append((n, n2))
        for n, n2 in pairs:
            p_map = {}
            for i, lp in enumerate(n.get_right_plugs()):
                network.remove_edge(lp)
                p_map[i] = lp
            for lp in n2.get_left_plugs():
                network.remove_edge(lp)
            for i, p in enumerate(n2.get_right_plugs()):
                p: Plug = p
                rp = p.edge.right_plug
                network.remove_edge(p)
                network.add_edge(p_map[i], rp)
            network.remove_simple(n2)
            network.change_node(n, self.type_result)


class HaarIntegration(TNComputation):
    def __init__(self, n_steps):
        self.two_haar = TwoHaarIntegration()
        self.four_haar = FourHaarIntegration()
        self.n_steps = n_steps

    def compute(self, networks: TensorNetworks) -> TensorNetworks:
        result = networks
        for s in range(self.n_steps):
            result = self.step(result)
        return result

    def step(self, networks: TensorNetworks):
        result = TensorNetworks()
        for i, n in enumerate(networks.networks):
            coeff = networks.coefficients[i]
            network: TensorNetwork = n
            self.integrate_one(result, coeff, network)
        return result

    def integrate_one(self, result: TensorNetworks, coeff: Coefficient, network: TensorNetwork):
        network.reduce()
        if len(network.group_map) == 0:
            result.add(coeff, network)
            return
        g_id, gates = sorted(network.group_map.items(), reverse=True)[0]
        if gates[0].type != Type.UNITARY:
            result.add(coeff, network)
            return
        if len(gates) == 2:
            coefficient = coeff.copy()
            factors, network = self.two_haar.integrate(network, g_id)
            coefficient.extend(factors)
            result.add(coefficient, network)
        elif len(gates) == 4:
            pairs = self.four_haar.integrate(network, g_id)
            for factors, network in pairs:
                coefficient = coeff.copy()
                coefficient.extend(factors)
                result.add(coefficient, network)
        return


class TwoHaarIntegration:
    def integrate(self, network: TensorNetwork, g_id):
        u, udagger = network.group_map[g_id]
        ys = GateUtil.get_ys(u)
        factors = []
        for y in ys:
            l = u.get_plug(Direction.Left, y)
            r = u.get_plug(Direction.Right, y)
            ld = udagger.get_plug(Direction.Left, y)
            rd = udagger.get_plug(Direction.Right, y)
            factor = self.do_integrate(network, l, r, ld, rd)
            factors.append(factor)
        network.remove_simple(u)
        network.remove_simple(udagger)
        return factors, network

    def do_integrate(self, network: TensorNetwork, l, r, ld, rd):
        lp = l.edge.left_plug
        rp = r.edge.right_plug
        ldp = ld.edge.left_plug
        rdp = rd.edge.right_plug
        network.remove_edge(l)
        network.remove_edge(r)
        network.remove_edge(ld)
        network.remove_edge(rd)
        if PlugUtil.connected(l, rd) and PlugUtil.connected(r, ld):
            factor = [Factor.D]
        elif PlugUtil.connected(l, rd):
            network.add_edge(lp, rp)
            factor = []
        elif PlugUtil.connected(r, ld):
            network.add_edge(lp, rdp)
            factor = []
        else:
            network.add_edge(ldp, rp)
            network.add_edge(lp, rdp)
            factor = [Factor.DF]
        return factor


class FourHaarIntegration:
    def __init__(self):
        self.integrators = [ParallelIntegrator(), ParallelIntegrator(True), CrossIntegrator(), CrossIntegrator(True)]
        self.post_processors = []

    def integrate(self, network: TensorNetwork, g_id):
        results = []
        hists = [Factor.HIST1, Factor.HIST2, Factor.HIST3, Factor.HIST4]
        for index, integrator in enumerate(self.integrators):
            net = network.copy()
            factor, net = integrator.integrate(net, g_id)
            factor.append(hists[index])
            results.append((factor, net))
            for post_process in self.post_processors:
                post_process.run(net)
        return results


class PostProcess(ABC):
    @abstractmethod
    def run(self, network):
        pass


class DrawNetwork(PostProcess):
    def run(self, network):
        network.draw()
        plt.show()


class Integrator(ABC):
    def __init__(self, swap=False):
        self.swap = swap

    def integrate(self, network: TensorNetwork, g_id):
        u1, udagger1, u2, udagger2 = network.group_map[g_id]
        ys1 = GateUtil.get_ys(u1)
        ys2 = GateUtil.get_ys(u2)
        factors = []
        for y1, y2 in zip(ys1, ys2):
            l1 = u1.get_plug(Direction.Left, y1)
            r1 = u1.get_plug(Direction.Right, y1)
            l1dagger = udagger1.get_plug(Direction.Left, y1)
            r1dagger = udagger1.get_plug(Direction.Right, y1)
            l2 = u2.get_plug(Direction.Left, y2)
            r2 = u2.get_plug(Direction.Right, y2)
            l2dagger = udagger2.get_plug(Direction.Left, y2)
            r2dagger = udagger2.get_plug(Direction.Right, y2)
            if self.swap:
                fs = self.do_integrate(network, l1, l2dagger, r1, r2dagger,
                                       l2, l1dagger, r2, r1dagger)
                factors.extend(fs)
            else:
                fs = self.do_integrate(network, l1, l1dagger, r1, r1dagger,
                                       l2, l2dagger, r2, r2dagger)
                factors.extend(fs)
        network.remove_simple(u1)
        network.remove_simple(u2)
        network.remove_simple(udagger1)
        network.remove_simple(udagger2)
        return factors, network

    def do_integrate(self, network: TensorNetwork,
                     l1: Plug, l1d: Plug,
                     r1: Plug, r1d: Plug,
                     l2: Plug, l2d: Plug,
                     r2: Plug, r2d: Plug):
        # pairs that becomes delta when integrated
        haar_pairs = self.get_haar_pairs(l1, l1d, r1, r1d, l2, l2d, r2, r2d)
        # pairs that connects after haar integration
        final_pairs, n_loop = PathUtil.find_outside_pairs(haar_pairs)
        map = {1: Factor.D, 2: Factor.D2, 3: Factor.D3, 4: Factor.D4}
        factors = self.initial_factors()
        if n_loop > 0:
            factors.append(map[n_loop])
        # remove edges in all plugs in integrated unitaries
        for p in haar_pairs:
            network.remove_edge(p[0])
            network.remove_edge(p[1])
        # connect edge between outside pairs
        for p in final_pairs:
            network.add_edge(p[0], p[1])
        return factors

    @abstractmethod
    def get_haar_pairs(self, l1: Plug, l1d: Plug, r1: Plug, r1d: Plug, l2: Plug, l2d: Plug, r2: Plug, r2d: Plug):
        return []

    @abstractmethod
    def initial_factors(self):
        return []


class PathUtil:
    @classmethod
    def find_outside_pairs(cls, haar_pairs):
        # path that is already connected
        pairs = cls.find_pairs(haar_pairs)
        left_plugs = []
        for p in haar_pairs:
            if not cls._included(p[0], pairs):
                left_plugs.append(p[0])
        final_pairs = []
        paths = []
        for left_p in left_plugs:
            left_p: Plug = left_p
            path = cls._add_path([left_p], left_p, pairs.copy(), haar_pairs.copy())
            paths.append(path)
            final_pairs.append((left_p.edge.left_plug, path[len(path) - 1].edge.right_plug))
        n_loop = 0
        for pair in pairs:
            path = cls._add_path([pair[0]], pair[0], pairs.copy(), haar_pairs.copy())
            if cls.contains(paths, path):
                continue
            paths.append(path)
            n_loop = n_loop + 1
        return final_pairs, n_loop

    @classmethod
    def contains(cls, paths, path):
        for path2 in paths:
            path_set = set(path2)
            for p in path:
                if p in path_set:
                    return True
        return False

    @classmethod
    def _add_path(cls, path, plug: Plug, pairs, haar_pairs):
        for i, p in enumerate(haar_pairs):
            if p[0] == plug:
                path.append(p[1])
                haar_pairs.pop(i)
                return cls._add_path(path, p[1], pairs, haar_pairs)
            elif p[1] == plug:
                path.append(p[0])
                haar_pairs.pop(i)
                return cls._add_path(path, p[0], pairs, haar_pairs)
        for i, p in enumerate(pairs):
            if p[0] == plug:
                if path.__contains__(p[1]):
                    return path
                pairs.pop(i)
                path.append(p[1])
                return cls._add_path(path, p[1], pairs, haar_pairs)
            if p[1] == plug:
                if path.__contains__(p[0]):
                    return path
                pairs.pop(i)
                path.append(p[0])
                return cls._add_path(path, p[0], pairs, haar_pairs)
        return path

    @classmethod
    def _included(cls, plug: Plug, pairs):
        for pair in pairs:
            if plug == pair[0] or plug == pair[1]:
                return True
        return False

    @classmethod
    def find_pairs(cls, haar_pairs):
        pairs = []
        for p in haar_pairs:
            for p2 in haar_pairs:
                if PlugUtil.connected(p2[1], p[0]):
                    pairs.append((p[0], p2[1]))
        return pairs


class ParallelIntegrator(Integrator):
    def get_haar_pairs(self, l1: Plug, l1d: Plug, r1: Plug, r1d: Plug, l2: Plug, l2d: Plug, r2: Plug, r2d: Plug):
        return [(l1, r1d), (l1d, r1), (l2, r2d), (l2d, r2)]

    def initial_factors(self):
        return [Factor.G]


class CrossIntegrator(Integrator):
    def get_haar_pairs(self, l1: Plug, l1d: Plug, r1: Plug, r1d: Plug, l2: Plug, l2d: Plug, r2: Plug, r2d: Plug):
        return [(l1, r1d), (l1d, r2), (l2, r2d), (l2d, r1)]

    def initial_factors(self):
        return [Factor.G, Factor.DF, Factor.MI]
