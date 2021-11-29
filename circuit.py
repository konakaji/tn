from core import *
from abc import ABC, abstractmethod


class LayerAppender(ABC):
    @abstractmethod
    def is_target(self, l_id):
        pass

    @abstractmethod
    def add(self, result, l_id, n_blk):
        pass


class TensorEncoderAppender(LayerAppender):
    def __init__(self, targets):
        self.targets = targets

    def is_target(self, l_id):
        return l_id in self.targets

    def add(self, result, l_id, n_blk):
        for j in range(2 * n_blk):
            loc = Location(l_id + 1, j, j)
            result.add_gate(Gate(loc, loc.id(), Type.ENC))


class ALTGenerator:
    @classmethod
    def generate(cls, n_blk, depth, obs_stt, obs_end, layer_appender: LayerAppender = None, revert=False):
        result = Circuit(2 * n_blk)
        result.add_gate(Gate(Location(0, 0, 2 * n_blk - 1),
                             Location(0, 0, 2 * n_blk - 1).id(), Type.INITIAL))
        sign = 0
        if revert:
            offset = 1
        for l in range(depth):
            if layer_appender is not None and layer_appender.is_target(l):
                layer_appender.add(result, l, n_blk)
                continue
            if (l + sign) % 2 == 0:
                cls.add_even_layer(result, l, n_blk)
            elif (l + sign) % 2 == 1:
                cls.add_odd_layer(result, l, n_blk)
        loc = Location(depth + 1, obs_stt, obs_end)
        result.add_observable(
            Gate(loc, "", Type.OBSERVABLE))
        return result

    @classmethod
    def add_even_layer(cls, result, l_id, n_blk):
        for j in range(n_blk):
            loc = Location(l_id + 1, 2 * j, 2 * j + 1)
            result.add_gate(
                Gate(loc, loc.id(), Type.UNITARY))

    @classmethod
    def add_odd_layer(cls, result, l_id, n_blk):
        loc = Location(l_id + 1, 0, 0)
        result.add_gate(Gate(loc, loc.id(), Type.UNITARY))
        for j in range(n_blk - 1):
            loc = Location(l_id + 1, 2 * j + 1, 2 * j + 2)
            result.add_gate(
                Gate(loc, loc.id(), Type.UNITARY))
        loc = Location(l_id + 1, 2 * n_blk - 1, 2 * n_blk - 1)
        result.add_gate(Gate(loc, loc.id(), Type.UNITARY))


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
        result.add(Coefficient(-1, []), net1)
        net2 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id)
        self._do_add_grad_avg(net2, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(Coefficient(2, []), net2)
        net3 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net3, grad_id, dagger=True)
        self._do_add_grad_avg(net3, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(Coefficient(-1, []), net3)
        for network in result.networks:
            network.transpile()
        return result

    def to_grad_avg(self, grad_id, mhalf) -> TensorNetworks:
        result = TensorNetworks()
        net1 = TensorNetwork(mhalf, self.b_height, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net1, grad_id)
        result.add(Coefficient(-1j, []), net1)
        net2 = TensorNetwork(mhalf, self.b_height, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id, dagger=True)
        result.add(Coefficient(1j, []), net2)
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
