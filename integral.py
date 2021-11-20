from core import *


class HaarIntegration:
    def __init__(self, n_steps):
        self.two_haar = TwoHaarIntegration()
        self.four_haar = FourHaarIntegration()
        self.n_steps = n_steps

    def integrate(self, networks: TensorNetworks) -> TensorNetworks:
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
            factors, network = self.two_haar.integrate(network, gates[0], gates[1])
            coefficient.extend(factors)
            result.add(coefficient, network)
        elif len(gates) == 4:
            pairs = self.four_haar.integrate(network, gates[0], gates[2], gates[1], gates[3])
        return


class TwoHaarIntegration:
    def integrate(self, network: TensorNetwork, u: Gate, udagger: Gate):
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
    def integrate(self, network: TensorNetwork, u1: Gate, u2: Gate, udagger1: Gate, udagger2: Gate):
        return []
