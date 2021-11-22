import networkx as nx
import math
import enum
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random


class Edge:
    def __init__(self, left_plug, right_plug):
        if left_plug.direction == Direction.Left \
                or right_plug.direction == Direction.Right:
            raise InvalidVariableException("The directions of plugs are invalid.")
        left_plug.edge = self
        right_plug.edge = self
        self.left_plug: Plug = left_plug
        self.right_plug: Plug = right_plug
        self.id = self.__hash__()

    def __hash__(self) -> int:
        if self.left_plug is None or self.right_plug is None:
            return 0
        return 31 + 13 * self.left_plug.__hash__() \
               + 17 * self.right_plug.__hash__()

    def __eq__(self, o: object) -> bool:
        if type(o) is not Edge:
            return False
        return self.__hash__() == o.__hash__()

    def detach(self):
        self.left_plug.edge = None
        self.right_plug.edge = None


class Factor(enum.Enum):
    # D = 2^m
    D = ("D", 1, 1)
    D2 = ("D^2", 2, 2)
    D3 = ("D^3", 3, 3)
    D4 = ("D^4", 4, 4)
    DF = ("1/D", 5, -1)
    MI = ("-", 6, 0)
    G = ("1/D^2-1", 7, 0)
    HIST1 = ("H1", 1, 0)
    HIST2 = ("H2", 2, 0)
    HIST3 = ("H3", 3, 0)
    HIST4 = ("H4", 4, 0)

    def __init__(self, label, id, d_count):
        self.label = label
        self.id = id
        self.d_count = d_count


class Direction(enum.Enum):
    Left = ("left", 23)
    Right = ("right", 31)

    def __init__(self, n, hash):
        self.n = n
        self.hash = hash

    def invert(self):
        if self == Direction.Left:
            return Direction.Right
        return Direction.Left

    def __hash__(self):
        return self.hash


class Coefficient:
    def __init__(self, digit, factors):
        self.factors = factors
        self.digit = digit

    def add(self, factor):
        self.factors.append(factor)

    def multiply(self, v):
        self.digit = self.digit * v

    def extend(self, factors):
        self.factors.extend(factors)

    def copy(self):
        result = Coefficient(self.digit, [])
        for f in self.factors:
            result.add(f)
        return result


class Plug:
    def __init__(self, j, direction, node_id):
        self.j = j
        self.direction = direction
        self.edge = None
        self.node_id = node_id
        self.id = random.randint(0, 10000000)

    def __hash__(self) -> int:
        return 13 * (self.j + 1) + 17 * self.node_id + self.direction.__hash__()

    def __eq__(self, o: object) -> bool:
        if type(o) is not Plug:
            return False
        return self.__hash__() == o.__hash__()


class Type(enum.Enum):
    UNITARY = ("U", 1)
    UNINTEGRABLE_UNITARY = ("W_", 2)
    OBSERVABLE = ("O", 3)
    GRAD = ("W", 4)  # the node where the gradient is computed
    UrWUr = ("UrWUr", 5)
    INITIAL = ("|0>", 6)

    def __init__(self, n, hash):
        self.n = n
        self.hash = hash

    def __hash__(self):
        return self.hash


class Location:
    def __init__(self, x, y_start, y_end):
        self.x = x
        self.y_start = y_start
        self.y_end = y_end

    def copy(self):
        return Location(self.x, self.y_start, self.y_end)

    def id(self):
        return self.x + 1 / (self.y_start + 1)

    def __hash__(self) -> int:
        return 13 * self.x + 17 * self.y_start + 31 * self.y_end

    def __eq__(self, o: object) -> bool:
        if type(o) is not Location:
            return False
        return self.__hash__() == o.__hash__()

    def __repr__(self) -> str:
        return "({}, {}-{})".format(self.x, self.y_start, self.y_end)


class GateUtil:
    @classmethod
    def get_ys(cls, gate):
        return [y for y in range(gate.get_location().y_start, gate.get_location().y_end + 1)]


class InvalidVariableException(Exception):
    pass


class PlugUtil:
    @classmethod
    def connected(cls, plug1, plug2):
        if plug1.direction == Direction.Left or plug2.direction == Direction.Right:
            raise InvalidVariableException("")
        if plug1.edge is None:
            return False
        if plug1.edge.right_plug == plug2:
            return True
        return False


class Gate:
    def __init__(self, location: Location, group_id, t: Type, dagger=False):
        self._location = location
        self.group_id = group_id
        self.type = t
        self.dagger = dagger
        self.id = self.__hash__()
        self.plugs = {Direction.Left: self._left_plugs(),
                      Direction.Right: self._right_plugs()}

    def copy_to(self, loc):
        return Gate(loc, self.group_id, self.type, dagger=self.dagger)

    def copy(self):
        return Gate(self.get_location().copy(), self.group_id, self.type, self.dagger)

    def get_location(self):
        return self._location.copy()

    def get_left_plugs(self):
        return self.plugs.get(Direction.Left)

    def get_right_plugs(self):
        return self.plugs.get(Direction.Right)

    def get_plug(self, direction: Direction, j):
        for plug in self.plugs[direction]:
            if plug.j == j:
                return plug
        return None

    def conjugate(self):
        c = self.copy_to(self._location.copy())
        c.dagger = not c.dagger
        return c

    def get_connectable(self, plug):
        for p in self.plugs[plug.direction.invert()]:
            if p.j == plug.j:
                return p
        return None

    def _left_plugs(self):
        if self.type == Type.INITIAL and not self.dagger:
            return []
        return [Plug(j, Direction.Left, self.id) for j in range(self._location.y_start, self._location.y_end + 1)]

    def _right_plugs(self):
        if self.type == Type.INITIAL and self.dagger:
            return []
        return [Plug(j, Direction.Right, self.id) for j in range(self._location.y_start, self._location.y_end + 1)]

    def __hash__(self) -> int:
        dag = 0
        if self.dagger:
            dag = 1
        return 13 * self._location.__hash__() + 17 * self.group_id.__hash__() \
               + 37 * dag + 41 * self.type.__hash__()

    def __eq__(self, o: object) -> bool:
        if type(o) is not Gate:
            return False
        return self.__hash__() == o.__hash__()

    def __repr__(self) -> str:
        dagger = ""
        if self.dagger:
            dagger = "†"
        return "{}{}".format(self.type.n, dagger)


class TensorNetwork:
    def __init__(self, mhalf, b_height, depth):
        self.mhalf = mhalf
        self.b_height = b_height
        self.depth = depth
        self.node_map = {}
        self.group_map = {}
        self.edge_map = {}

    def copy(self):
        result = TensorNetwork(self.mhalf, self.b_height, self.depth)
        for g in self.node_map.values():
            result.add_node(g.copy())
        for e in self.edge_map.values():
            edge: Edge = e
            lp_: Plug = edge.left_plug
            rp_: Plug = edge.right_plug
            lg: Gate = result.node_map[lp_.node_id]
            rg: Gate = result.node_map[rp_.node_id]
            lp, rp = None, None
            for p in lg.get_right_plugs():
                if p.j == lp_.j:
                    lp = p
                    break
            for p in rg.get_left_plugs():
                if p.j == rp_.j:
                    rp = p
            result.add_edge(lp, rp)
        return result

    def change_node(self, node: Gate, t: Type):
        self.node_map.pop(node.id)
        node.type = t
        node.id = node.__hash__()
        left_plugs = node.get_left_plugs()
        right_plugs = node.get_right_plugs()
        for p, p2 in zip(left_plugs, right_plugs):
            p.node_id = node.id
            p2.node_id = node.id
        self.node_map[node.id] = node

    def add_node(self, gate: Gate):
        self.node_map[gate.id] = gate
        if gate.type != Type.UNITARY:
            return
        if gate.group_id not in self.group_map:
            self.group_map[gate.group_id] = []
        self.group_map[gate.group_id].append(gate)

    def add_edge(self, lp: Plug, rp: Plug):
        edge = Edge(lp, rp)
        self.edge_map[edge.id] = edge

    def nodes(self):
        return self.node_map.values()

    def reduce(self):
        for node in self.nodes():
            if node.type != Type.UNITARY or len(node.get_right_plugs()) == 0:
                continue
            n_id = None
            add = True
            for p in node.get_right_plugs():
                if n_id is None:
                    n_id = p.edge.right_plug.node_id
                elif n_id != p.edge.right_plug.node_id:
                    add = False
                    break
                if node.group_id != self.node_map[n_id].group_id:
                    add = False
                    break
            if add:
                self.remove(node, self.node_map[n_id])
                self.reduce()
                return
        return

    def remove_edge(self, plug: Plug):
        if plug.edge is not None:
            self.edge_map.pop(plug.edge.id)
            plug.edge.detach()

    def remove_simple(self, node: Gate):
        self.node_map.pop(node.id)
        members = []
        if node.group_id not in self.group_map:
            return
        for n in self.group_map[node.group_id]:
            if n.id == node.id:
                continue
            members.append(n)
        if len(members) == 0:
            self.group_map.pop(node.group_id)
        else:
            self.group_map[node.group_id] = members

    def remove(self, l_node: Gate, r_node: Gate):
        left_map = {}
        right_map = {}
        for i, p in enumerate(l_node.get_left_plugs()):
            plug: Plug = p
            if plug.edge is None:
                continue
            lp: Plug = plug.edge.left_plug
            self.remove_edge(plug)
            left_map[i] = lp
        for p in l_node.get_right_plugs():
            self.remove_edge(p)
        for i, p in enumerate(r_node.get_right_plugs()):
            plug: Plug = p
            if plug.edge is None:
                continue
            rp: Plug = plug.edge.right_plug
            self.remove_edge(plug)
            right_map[i] = rp
        for j, lp in left_map.items():
            rp = right_map[j]
            self.add_edge(lp, rp)
        self.node_map.pop(l_node.id)
        self.node_map.pop(r_node.id)
        members = []
        group_id = l_node.group_id
        for n in self.group_map[group_id]:
            if n.id == l_node.id or n.id == r_node.id:
                continue
            members.append(n)
        if len(members) == 0:
            self.group_map.pop(group_id)
        else:
            self.group_map[group_id] = members

    def transpile(self):
        nodes = sorted(self.nodes(), key=lambda n: n.get_location().x)
        for n in nodes:
            for plug in n.get_right_plugs():
                for n2 in nodes:
                    if n.get_location().x >= n2.get_location().x:
                        continue
                    right = n2.get_connectable(plug)
                    if right is not None:
                        edge = Edge(plug, right)
                        self.edge_map[edge.id] = edge
                        break

    def draw(self, grid_width=1, space=0.3, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        ax.set_xlim([-0.1, self.depth])
        ax.set_ylim([0, self.b_height + 0.5])
        netx = nx.DiGraph()
        right_plugs = []
        for node in self.nodes():
            for plug in node.get_left_plugs():
                netx.add_node(plug.id, pos=(node.get_location().x, plug.j + 0.5))
            for plug in node.get_right_plugs():
                right_plugs.append(plug)
                netx.add_node(plug.id,
                              pos=(node.get_location().x + grid_width - space, plug.j + 0.5 * grid_width))
            rect = patches.Rectangle(
                xy=(grid_width * node.get_location().x, space + grid_width * node.get_location().y_start),
                width=grid_width - space,
                height=grid_width * (node.get_location().y_end - node.get_location().y_start + 0.9),
                linewidth=1, edgecolor='black', facecolor='none')
            ax.add_patch(rect)
            plt.text(node.get_location().x + 0.1, node.get_location().y_start + 0.5, str(node))
        for plug in right_plugs:
            if plug.edge is not None:
                netx.add_edge(plug.id, plug.edge.right_plug.id)
        pos = nx.get_node_attributes(netx, 'pos')
        nx.draw_networkx(netx, pos, node_size=10,
                         with_labels=False, connectionstyle="arc3,rad=0.1")

    def __eq__(self, o: object) -> bool:
        if type(o) is not TensorNetwork:
            return False
        return self.__hash__() == o.__hash__()

    def __hash__(self):
        result = 0
        for i, item in enumerate(sorted(self.node_map.items())):
            k, gate = item
            result = result + gate.__hash__()
        for i, item in enumerate(sorted(self.edge_map.items())):
            k, edge = item
            result = result + edge.__hash__() * (i + 1)
        return result


class TensorNetworks:
    def __init__(self):
        self.coefficients = []
        self.networks = []

    def add(self, coeff, network):
        self.coefficients.append(coeff)
        self.networks.append(network)

    def draw(self, figsize=(10, 10)):
        fig = plt.figure(figsize=figsize)
        length = int(math.sqrt(len(self.coefficients))) + 1
        for i, coeff in enumerate(self.coefficients):
            ax = fig.add_subplot(length, length, i + 1)
            self.networks[i].draw(ax=ax)
