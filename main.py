from circuit import *
from core import *
from integral import *

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
    integral = HaarIntegration(2)
    result = integral.integrate(tn)
    result.draw(figsize=(10, 10))
    plt.show()