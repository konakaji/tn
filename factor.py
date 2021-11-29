from core import *


class CoefficientInterpreter:
    @classmethod
    def interpret(cls, coefficient):
        hist_map = {Factor.HIST1: 1, Factor.HIST2: 2, Factor.HIST3: 3, Factor.HIST4: 4}
        histories = []
        for f in coefficient.factors:
            if f in hist_map:
                histories.append(hist_map[f])
        factors = sorted(coefficient.factors, key=lambda f: f.id)
        d_count = 0
        g_count = 0
        digit = coefficient.digit
        for f in factors:
            f: Factor = f
            d_count = d_count + f.d_count
            if f == Factor.G:
                g_count = g_count + 1
            elif f == Factor.G2:
                g_count = g_count + 2
            elif f == Factor.MI:
                digit = digit * -1
        return DecayFactor(digit, d_count, g_count, histories)


class DecayFactor:
    def __init__(self, digit, d_count, g_count, histories):
        self.digit = digit
        self.d_count = d_count
        self.g_count = g_count
        self.histories = histories

    def approximate_d_count(self):
        if self.digit == 0:
            return 0
        return self.d_count - 2 * self.g_count

    def is_appendable(self, d):
        d: DecayFactor = d
        return self.d_count == d.d_count and self.g_count == d.g_count

    def append(self, d):
        if not self.is_appendable(d):
            raise InvalidVariableException("not appendable")
        self.histories = []
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


class FactorUtil:
    @classmethod
    def get_max_dcount(cls, decay_factors):
        result = -1000000
        for df in decay_factors:
            df: DecayFactor = df
            ac = df.approximate_d_count()
            if ac != 0 and ac > result:
                result = ac
        return result
