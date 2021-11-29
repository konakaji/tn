from tn.core import *
import os, matplotlib.pyplot as plt, shutil


class CoefficientInterpreter:
    @classmethod
    def interpret(cls, coefficient):
        hist_map = {Factor.HIST1: 1, Factor.HIST2: 2, Factor.HIST3: 3, Factor.HIST4: 4,
                    Factor.HIST_A: "a", Factor.HIST_B: "b", Factor.HIST_C: "c"}
        histories = []
        for f in coefficient.factors:
            if f in hist_map:
                histories.append(f.label)
        factors = sorted(coefficient.factors, key=lambda f: f.id)
        d_count = 0
        g_count = 0
        g2_count = 0
        digit = coefficient.digit
        for f in factors:
            f: Factor = f
            d_count = d_count + f.d_count
            if f == Factor.G:
                g_count = g_count + 1
            elif f == Factor.G2:
                g2_count = g2_count + 1
            elif f == Factor.MI:
                digit = digit * -1
        return CoefficientReduced(digit, d_count, g_count, g2_count, histories)


class CoefficientReduced:
    def __init__(self, digit, d_count, g_count, g2_count, histories):
        self.digit = digit
        self.d_count = d_count
        self.g_count = g_count
        self.g2_count = g2_count
        self.histories = histories

    def approximate_d_count(self):
        if self.digit == 0:
            return 0
        return self.d_count - 2 * self.g_count - 4 * self.g2_count

    def is_appendable(self, d):
        d: CoefficientReduced = d
        return self.d_count == d.d_count and self.g_count == d.g_count

    def append(self, d):
        if not self.is_appendable(d):
            raise InvalidVariableException("not appendable")
        self.histories = []
        self.digit = self.digit + d.digit

    def copy(self):
        return CoefficientReduced(self.digit, self.d_count, self.g_count, self.g2_count, self.histories)

    def __repr__(self):
        return "{} x 2^({}m)/(2^(2m)-1)^{} (2^(4m)-1)^{}".format(self.digit, self.d_count,
                                                                 self.g_count, self.g2_count)


class CoefficientMerger:
    @classmethod
    def merge(cls, coefficients):
        excludes = set()
        results = []
        coefficients = cls.copy(coefficients)
        for c in coefficients:
            c: CoefficientReduced = c
            if c in excludes:
                continue
            excludes.add(c)
            cls.do_merge(c, coefficients, excludes)
            results.append(c)
        return results

    @classmethod
    def copy(cls, coefficients):
        results = []
        for c in coefficients:
            results.append(c.copy())
        return results

    @classmethod
    def do_merge(cls, c, coefficients, excludes):
        for c2 in coefficients:
            if c2 not in excludes and c.is_appendable(c2):
                excludes.add(c2)
                c.append(c2)


class CoefficientUtil:
    @classmethod
    def get_max_dcount(cls, decay_factors):
        result = -1000000
        for df in decay_factors:
            df: CoefficientReduced = df
            ac = df.approximate_d_count()
            if ac != 0 and ac > result:
                result = ac
        return result


class ReportBuilder:
    def __init__(self, result: TensorNetworks):
        self.result = result
        self.network_map = None
        self.merged_map = None
        self.d_map = None

    def build(self):
        self.network_map = self._build_network_map()
        self.merged_map = self._merge_coefficient()
        self.d_map = self._create_max_d_count()
        return self

    def to_html(self, path):
        image_path = "{}/images".format(path)
        detail_path = "{}/details".format(path)
        if not os.path.exists(path):
            os.mkdir(path)
            os.mkdir(image_path)
            os.mkdir(detail_path)
        else:
            shutil.rmtree(path)
            os.mkdir(path)
            os.mkdir(image_path)
            os.mkdir(detail_path)
        for network, coeffs in self.network_map.items():
            plt.clf()
            network.draw()
            plt.savefig("{}/{}".format(image_path, network.__hash__()))
        for network in self.result.networks:
            with open("{}/{}.html".format(detail_path, network.__hash__()), "w") as f:
                f.write("<html>")
                f.write("<head>")
                f.write(
                    '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta1/'
                    'dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-g'
                    'iJF6kkoqNQ00vy+HMDP7azOuL0xtbfIcaT9wjKHr8RbDVddVHyTfAAsrekwKmP1"'
                    ' crossorigin="anonymous">')
                f.write("</head>")
                f.write("<h1>Detail</h1>\n")
                f.write("<img src='../../{}/{}.png' width=400></img>\n"
                        .format(image_path, network.__hash__()))
                f.write("<table class='table table-bordered'>")
                coeffs = self.network_map[network]
                f.write("<tr><th>Coefficient</th><th>approximate_d_count</th><th>Histories</th></tr>\n")
                for coeff in coeffs:
                    f.write("<tr><td>{}</td><td>{}</td><td>{}</td></tr>\n".format(coeff, coeff.approximate_d_count(),
                                                                                  coeff.histories))
                f.write("</table>\n")
                f.write("</html>")
        with open("{}/index.html".format(path), "w") as f:
            f.write("<html>")
            f.write("<head>")
            f.write(
                '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.0-beta1/'
                'dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-g'
                'iJF6kkoqNQ00vy+HMDP7azOuL0xtbfIcaT9wjKHr8RbDVddVHyTfAAsrekwKmP1"'
                ' crossorigin="anonymous">')
            f.write("</head>")
            f.write("<table class='table table-bordered'>")
            f.write("<tr><th>Network</th><th>d_count</th><th>action</th></tr>")
            for network, d_count in sorted(self.d_map.items(), key=lambda v: -v[1]):
                f.write("<tr>")
                f.write("<td><img src='../{}/{}.png' width=300></img></td>\n"
                        .format(image_path, network.__hash__()))
                f.write("<td>{}</td>".format(d_count))
                f.write("<td><a href='details/{}.html'>detail</a></td>".format(network.__hash__()))
                f.write("<tr>")
            f.write("</table>")
            f.write("</html>")

    def _build_network_map(self):
        network_map = {}
        for i, network in enumerate(self.result.networks):
            coeff = self.result.coefficients[i]
            if network not in network_map:
                network_map[network] = []
            network_map[network].append(CoefficientInterpreter.interpret(coeff))
        return network_map

    def _merge_coefficient(self):
        result = {}
        for network, coeffs in self.network_map.items():
            fs = []
            for coeff in coeffs:
                fs.append(coeff)
            result[network] = CoefficientMerger.merge(fs)
        return result

    def _create_max_d_count(self):
        result = {}
        for network, merged in self.merged_map.items():
            result[network] = CoefficientUtil.get_max_dcount(merged)
        return result
