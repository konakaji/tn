import networkx as nx
import math
import enum
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random


class Node:
    def __init__(self, n_qubit, id):
        self.id = id
        self.nqubit = n_qubit


class Edge:
    def __init__(self, left_plug, right_plug):
        left_plug.edge = self
        right_plug.edge = self
        self.left_plug = left_plug
        self.right_plug = right_plug
        self.id = random.randint(0, 100000000)

    def detach(self):
        self.left_plug.edge = None
        self.right_plug.edge = None


class Direction(enum.Enum):
    Left = "left"
    Right = "right"

    def invert(self):
        if self == Direction.Left:
            return Direction.Right
        return Direction.Left


class Plug:
    def __init__(self, j, direction, node_id):
        self.j = j
        self.direction = direction
        self.edge = None
        self.node_id = node_id
        self.id = random.randint(0, 10000000)


class Type(enum.Enum):
    UNITARY = "U"
    UNINTEGRABLE_UNITARY = "Uw"
    OBSERVABLE = "O"
    GRAD = "W"  # the node where the gradient is computed
    INITIAL = "|0>"


class Location:
    def __init__(self, x, y_start, y_end):
        self.x = x
        self.y_start = y_start
        self.y_end = y_end

    def copy(self):
        return Location(self.x, self.y_start, self.y_end)

    def id(self):
        return self.x + 1 / (self.y_start + 1)

    def __repr__(self) -> str:
        return "({}, {}-{})".format(self.x, self.y_start, self.y_end)


class Gate:
    def __init__(self, location: Location, group_id, t: Type, id=None, dagger=False):
        self._location = location
        self.group_id = group_id
        self.type = t
        if id is None:
            id = random.randint(0, 1000000000)
        self.id = id
        self.dagger = dagger
        self.plugs = {Direction.Left: self._left_plugs(),
                      Direction.Right: self._right_plugs()}

    def copy_to(self, loc):
        return Gate(loc, self.group_id, self.type, dagger=self.dagger)

    def copy(self):
        return Gate(self.get_location().copy(), self.group_id, self.type, self.id, self.dagger)

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

    def __repr__(self) -> str:
        dagger = ""
        if self.dagger:
            dagger = "†"
        return "{}{}".format(self.type.value, dagger)


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

    def remove(self, l_node: Gate, r_node: Gate):
        left_map = {}
        right_map = {}
        for p in l_node.get_left_plugs():
            plug: Plug = p
            lp: Plug = plug.edge.left_plug
            self.edge_map.pop(plug.edge.id)
            plug.edge.detach()
            left_map[lp.j] = lp
        for p in l_node.get_right_plugs():
            self.edge_map.pop(p.edge.id)
        for p in r_node.get_right_plugs():
            plug: Plug = p
            rp: Plug = plug.edge.right_plug
            self.edge_map.pop(plug.edge.id)
            plug.edge.detach()
            right_map[rp.j] = rp
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
        nx.draw_networkx(netx, pos, node_size=10, with_labels=False)


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


class Circuit:
    def __init__(self, b_height):
        self.gates = []
        self.b_height = b_height
        self.observable = None

    def add_gate(self, gate: Gate):
        self.gates.append(gate)

    def add_observable(self, observable: Gate):
        self.observable = observable

    def to_grad_var(self, grad_id, mhalf) -> TensorNetworks:
        result = TensorNetworks()
        net1 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net1, grad_id)
        self._do_add_grad_avg(net1, grad_id, y_offset=self.b_height + 1)
        result.add(-1, net1)
        net2 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id)
        self._do_add_grad_avg(net2, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(2, net2)
        net3 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net3, grad_id, dagger=True)
        self._do_add_grad_avg(net3, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(-1, net3)
        for network in result.networks:
            network.transpile()
        return result

    def to_grad_avg(self, grad_id, mhalf) -> TensorNetworks:
        result = TensorNetworks()
        net1 = TensorNetwork(mhalf, self.b_height, self.observable.location.x * 2 + 1)
        self._do_add_grad_avg(net1, grad_id)
        result.add(-1j, net1)
        net2 = TensorNetwork(mhalf, self.b_height, self.observable.location.x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id, dagger=True)
        result.add(1j, net2)
        for network in result.networks:
            network.transpile()
        return result

    def _do_add_grad_avg(self, result: TensorNetwork, grad_id, dagger=False, y_offset=0):
        for gate in self.gates:
            loc = Location(gate.get_location().x,
                           y_offset + gate.get_location().y_start,
                           y_offset + gate.get_location().y_end)
            if grad_id == gate.group_id:
                t = Type.GRAD
                if dagger:
                    t = Type.UNINTEGRABLE_UNITARY
                result.add_node(Gate(loc, gate.group_id, t))
            else:
                result.add_node(gate.copy_to(loc))
        o = self.observable.copy_to(Location(self.observable.get_location().x,
                                             self.observable.get_location().y_start + y_offset,
                                             self.observable.get_location().y_end + y_offset))
        result.add_node(o)
        offset = self.observable.get_location().x * 2
        for gate in reversed(self.gates):
            loc = Location(offset - gate.get_location().x,
                           y_offset + gate.get_location().y_start,
                           y_offset + gate.get_location().y_end)
            if grad_id == gate.group_id:
                t = Type.UNINTEGRABLE_UNITARY
                if dagger:
                    t = Type.GRAD
                result.add_node(Gate(loc, gate.group_id, t, dagger=True))
            else:
                g = gate.conjugate()
                g.location = loc
                result.add_node(g.copy_to(loc))
        return result


class HaarIntegration:
    def integrate(self, networks: TensorNetworks) -> TensorNetworks:
        for i, n in enumerate(networks.networks):
            coeff = networks.coefficients[i]
            network: TensorNetwork = n
            self.step(coeff, network)

    def step(self, coeff, network) -> TensorNetworks:
        network.reduce()
        for g_id, gates in sorted(network.group_map.items(), reverse=True):
            break


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
    tn = network.to_grad_var("a", 1)
    for net in tn.networks:
        net.reduce()
    tnn: TensorNetwork = tn.networks[0].copy()
    tnn.draw()
    # integral = HaarIntegration()
    # integral.integrate(tn)
    # tn.draw()
    plt.show()
    #
    # network = TensorNetwork()
    # network.add_node(Gate(Location(0, 0, 1), "hoge", Type.UNITARY, True))
    # network.add_node(Gate(Location(0, 2, 3), "fuga", Type.UNITARY, True))
    # network.add_node(Gate(Location(1, 0, 2), "fuga", Type.UNITARY, True))
    # network.add_node(Gate(Location(1, 3, 3), "fuga", Type.UNITARY, True))
    # network.draw()
    # plt.show()
# See PyCharm help at https://www.jetbrains.com/help/pycharm/
