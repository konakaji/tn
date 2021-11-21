from core import *


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
        result.add(Coefficient([-1]), net1)
        net2 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id)
        self._do_add_grad_avg(net2, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(Coefficient([2]), net2)
        net3 = TensorNetwork(mhalf, self.b_height * 2 + 1, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net3, grad_id, dagger=True)
        self._do_add_grad_avg(net3, grad_id, dagger=True, y_offset=self.b_height + 1)
        result.add(Coefficient([-1]), net3)
        for network in result.networks:
            network.transpile()
        return result

    def to_grad_avg(self, grad_id, mhalf) -> TensorNetworks:
        result = TensorNetworks()
        net1 = TensorNetwork(mhalf, self.b_height, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net1, grad_id)
        result.add(Coefficient([-1j]), net1)
        net2 = TensorNetwork(mhalf, self.b_height, self.observable.get_location().x * 2 + 1)
        self._do_add_grad_avg(net2, grad_id, dagger=True)
        result.add(Coefficient([1j]), net2)
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
