from circuit import *
from core import *
from computation import *


class CoefficientInterpreter:
    @classmethod
    def interpret(cls, coefficient):
        factors = sorted(coefficient.factors, key=lambda f: f.id)
        d_count = 0
        g_count = 0
        digit = coefficient.digit
        for f in factors:
            f: Factor = f
            d_count = d_count + f.d_count
            if f == Factor.G:
                g_count = g_count + 1
            elif f == Factor.MI:
                digit = digit * -1
        return DecayFactor(digit, d_count, g_count)


class DecayFactor:
    def __init__(self, digit, d_count, g_count):
        self.digit = digit
        self.d_count = d_count
        self.g_count = g_count

    def approximate_d_count(self):
        return self.d_count - 2 * self.g_count

    def is_appendable(self, d):
        d: DecayFactor = d
        return self.d_count == d.d_count and self.g_count == d.g_count

    def append(self, d):
        if not self.is_appendable(d):
            raise InvalidVariableException("not appendable")
        self.digit = self.digit + d.digit

    def __repr__(self):
        return "{} x 2^({}m)/(2^(2m)-1)^{}".format(self.digit, self.d_count, self.g_count)


class FactorMerger:
    @classmethod
    def merge(cls, factors):
        excludes = set()
        results = []
        for f in factors:
            f: DecayFactor = f
            if f in excludes:
                continue
            excludes.add(f)
            cls.do_merge(f, factors, excludes)
            results.append(f)
        return results

    @classmethod
    def do_merge(cls, f, factors, excludes):
        for f2 in factors:
            if f2 not in excludes and f.is_appendable(f2):
                excludes.add(f2)
                f.append(f2)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    network = Circuit(4)
    network.add_gate(Gate(Location(0, 0, 3), Location(0, 0, 3).id(), Type.INITIAL))
    network.add_gate(Gate(Location(1, 0, 1), Location(1, 0, 1).id(), Type.UNITARY))
    network.add_gate(Gate(Location(1, 2, 3), Location(1, 2, 3).id(), Type.UNITARY))
    network.add_gate(Gate(Location(2, 0, 0), Location(2, 0, 0).id(), Type.UNITARY))
    network.add_gate(Gate(Location(2, 1, 2), Location(2, 1, 2).id(), Type.UNITARY))
    network.add_gate(Gate(Location(2, 3, 3), Location(2, 3, 3).id(), Type.UNITARY))
    network.add_observable(Gate(Location(3, 1, 1), "", Type.OBSERVABLE))
    tn = network.to_grad_var(Location(1, 0, 1).id(), 1)
    computations = [HaarIntegration(2),
                    Multiply(Type.GRAD, Type.UNINTEGRABLE_UNITARY, False, True, Type.UrWUr, False),
                    Multiply(Type.UNINTEGRABLE_UNITARY, Type.GRAD, False, True, Type.UrWUr, False)]
    result = tn
    for computation in computations:
        result = computation.compute(result)
    network_map = {}
    for i, network in enumerate(result.networks):
        coeff = result.coefficients[i]
        if network not in network_map:
            network_map[network] = []
        network_map[network].append(coeff)
    ns = []
    for network, coeffs in network_map.items():
        fs = []
        for coeff in coeffs:
            f = CoefficientInterpreter.interpret(coeff)
            if f.approximate_d_count() == -2:
                ns.append(network)
            fs.append(f)
        print(network.__hash__(), FactorMerger.merge(fs))
